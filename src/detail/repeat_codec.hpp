#pragma once
#include "bytes.hpp"
#include "hash.hpp"
#include "predictor_v3.hpp"
#include "range.hpp"
#include "v4zip/dna.hpp"
#include <algorithm>
#include <array>
#include <cstdint>
#include <vector>

namespace v4zip::detail {

// Referential/Self-LZ + PredictorV3 hybrid for high-repeat / multi-strain / pan-genome DNA.
// Emits: [mode_tokens_zstd] + [residual_literals_range_coded]
class RepeatDnaCodec {
    static constexpr unsigned K = 16;
    static constexpr unsigned MIN_MATCH = 18;

public:
    static Bytes encode(const Bytes& dna, const DnaConfig& config) {
        const std::size_t n = dna.size();
        // Hash table: canonical 16-mer -> last 2 positions (packed (pos << 2) | transform)
        const std::size_t tbl_size = std::size_t{1} << 22; // 4M buckets = 64 MB
        struct Entry {
            std::uint64_t tag = UINT64_MAX;
            std::array<std::uint64_t, 2> pos{};
        };
        std::vector<Entry> table(tbl_size);

        Bytes tokens;      // 1-bit/byte flags + varints: 0=lit_run(len), 1=match(flags, dist, len)
        Bytes literals;
        DnaConfig lit_cfg;
        lit_cfg.format_version = config.format_version;

        const std::uint64_t mask = word_mask(K);
        std::uint64_t w = 0;
        std::size_t lit_start = 0;
        std::size_t i = 0;

        auto flush_lits = [&](std::size_t upto) {
            if (upto <= lit_start) return;
            const std::size_t run = upto - lit_start;
            tokens.push_back(0);
            varint(tokens, run);
            literals.insert(literals.end(), dna.begin() + static_cast<std::ptrdiff_t>(lit_start),
                            dna.begin() + static_cast<std::ptrdiff_t>(upto));
            lit_start = upto;
        };

        auto insert_kmer = [&](std::size_t end_pos, std::uint64_t word) {
            const auto canon = canonical_word(word, K);
            auto& row = table[hash_word(canon.word) & (tbl_size - 1)];
            if (row.tag != canon.word) { row.tag = canon.word; row.pos = {}; }
            row.pos[1] = row.pos[0];
            row.pos[0] = (std::uint64_t(end_pos) << 2) | canon.transform;
        };

        while (i < n) {
            w = ((w << 2) | (dna[i] & 3u)) & mask;
            if (i + 1 - lit_start < K) { ++i; continue; }

            const auto canon = canonical_word(w, K);
            const auto& row = table[hash_word(canon.word) & (tbl_size - 1)];
            std::size_t best_len = 0, best_src = 0;
            unsigned best_flags = 0; // bit0: rev_dir, bit1: comp

            if (row.tag == canon.word) {
                for (auto packed : row.pos) {
                    if (!packed) continue;
                    const std::size_t prev_end = static_cast<std::size_t>(packed >> 2);
                    const unsigned g = canon.transform ^ static_cast<unsigned>(packed & 3u);
                    const bool rev_dir = (g & 1u) != 0;
                    const unsigned comp = (g & 2u) ? 3u : 0u;
                    const std::size_t cur_start = i + 1 - K;

                    // Compute starting source index of the K-mer in history
                    if (!rev_dir) {
                        if (prev_end < K || prev_end - K >= cur_start) continue;
                        const std::size_t src_start = prev_end - K;
                        bool ok = true;
                        for (unsigned d = 0; d < K; ++d) {
                            if ((dna[src_start + d] ^ comp) != dna[cur_start + d]) { ok = false; break; }
                        }
                        if (!ok) continue;
                        std::size_t mlen = K;
                        while (cur_start + mlen < n && src_start + mlen < cur_start &&
                               (dna[src_start + mlen] ^ comp) == dna[cur_start + mlen]) {
                            ++mlen;
                        }
                        if (mlen > best_len) {
                            best_len = mlen;
                            best_src = src_start;
                            best_flags = comp ? 2u : 0u;
                        }
                    } else {
                        if (prev_end < K || prev_end - 1 >= cur_start) continue;
                        const std::size_t src_start = prev_end - 1;
                        bool ok = true;
                        for (unsigned d = 0; d < K; ++d) {
                            if ((dna[src_start - d] ^ comp) != dna[cur_start + d]) { ok = false; break; }
                        }
                        if (!ok) continue;
                        std::size_t mlen = K;
                        while (cur_start + mlen < n && src_start >= mlen &&
                               (dna[src_start - mlen] ^ comp) == dna[cur_start + mlen]) {
                            ++mlen;
                        }
                        if (mlen > best_len) {
                            best_len = mlen;
                            best_src = src_start;
                            best_flags = 1u | (comp ? 2u : 0u);
                        }
                    }
                }
            }

            if (best_len >= MIN_MATCH) {
                const std::size_t cur_start = i + 1 - K;
                flush_lits(cur_start);
                tokens.push_back(static_cast<std::uint8_t>(1u + best_flags));
                varint(tokens, cur_start - best_src);
                varint(tokens, best_len);
                // Index k-mers along the skipped match (stride 4) so future regions can match it too
                for (std::size_t j = i; j < cur_start + best_len; ++j) {
                    if (j > i) w = ((w << 2) | (dna[j] & 3u)) & mask;
                    if (((j + 1) & 3u) == 0 || j == i) insert_kmer(j + 1, w);
                }
                i = cur_start + best_len;
                lit_start = i;
                continue;
            }

            insert_kmer(i + 1, w);
            ++i;
        }
        flush_lits(n);

        // Compress literals with PredictorV3 and tokens with Zstd-19
        PredictorV3 model(DnaStrategy::coding_large);
        RangeEncoder enc;
        for (std::size_t j = 0; j < literals.size(); ++j) {
            const auto cdf = model.predict(literals, j);
            enc.put(literals[j], cdf);
            model.observe(literals[j], j);
        }
        const Bytes enc_lits = enc.finish();
        const Bytes enc_toks = pack(tokens);

        Bytes out;
        varint(out, literals.size());
        varint(out, tokens.size());
        varint(out, enc_toks.size());
        out.insert(out.end(), enc_toks.begin(), enc_toks.end());
        out.insert(out.end(), enc_lits.begin(), enc_lits.end());
        return out;
    }

    static Bytes decode(const Bytes& payload, std::size_t expected_dna) {
        std::size_t pos = 0;
        const std::size_t lit_count = read_size(payload, pos);
        const std::size_t tok_raw_size = read_size(payload, pos);
        const std::size_t tok_enc_size = read_size(payload, pos);
        if (tok_enc_size > payload.size() - pos) fail("truncated repeat tokens");
        const Bytes tokens = unpack(payload.data() + pos, tok_enc_size, tok_raw_size);
        pos += tok_enc_size;
        const Bytes enc_lits(payload.begin() + static_cast<std::ptrdiff_t>(pos), payload.end());

        if (lit_count > expected_dna) fail("invalid repeat literal count");
        PredictorV3 model(DnaStrategy::coding_large);
        RangeDecoder dec(enc_lits);
        Bytes literals;
        literals.reserve(lit_count);
        for (std::size_t j = 0; j < lit_count; ++j) {
            const auto base = dec.get(model.predict(literals, j));
            literals.push_back(static_cast<std::uint8_t>(base));
            model.observe(base, j);
        }
        dec.finish();

        Bytes out;
        out.reserve(expected_dna);
        std::size_t tp = 0, lp = 0;
        while (tp < tokens.size()) {
            const unsigned kind = tokens[tp++];
            if (kind == 0) {
                const std::size_t run = read_varint(tokens, tp);
                if (run > lit_count - lp || run > expected_dna - out.size()) fail("invalid repeat literal run");
                out.insert(out.end(), literals.begin() + static_cast<std::ptrdiff_t>(lp),
                           literals.begin() + static_cast<std::ptrdiff_t>(lp + run));
                lp += run;
            } else if (kind >= 1 && kind <= 4) {
                const unsigned flags = kind - 1;
                const bool rev_dir = (flags & 1u) != 0;
                const unsigned comp = (flags & 2u) ? 3u : 0u;
                const std::size_t dist = read_varint(tokens, tp);
                const std::size_t mlen = read_varint(tokens, tp);
                if (dist == 0 || dist > out.size() || mlen > expected_dna - out.size()) fail("invalid repeat match");
                const std::size_t src = out.size() - dist;
                if (!rev_dir) {
                    for (std::size_t k = 0; k < mlen; ++k) out.push_back(static_cast<std::uint8_t>(out[src + k] ^ comp));
                } else {
                    if (mlen > src + 1) fail("invalid reverse repeat match");
                    for (std::size_t k = 0; k < mlen; ++k) out.push_back(static_cast<std::uint8_t>(out[src - k] ^ comp));
                }
            } else {
                fail("invalid repeat token");
            }
        }
        if (lp != lit_count || out.size() != expected_dna) fail("repeat dna length mismatch");
        return out;
    }
};

} // namespace v4zip::detail
