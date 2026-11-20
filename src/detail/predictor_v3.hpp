#pragma once
#include "genome_profiler.hpp"
#include "hash.hpp"
#include "range.hpp"
#include "weights.hpp"
#include <algorithm>
#include <array>
#include <cstdint>
#include <vector>

namespace v4zip::detail {

struct CompactCounter {
    std::array<std::uint16_t, 4> c{};
    void observe(unsigned b, unsigned cap) {
        if (c[b] < 65000u) ++c[b];
        const unsigned sum = unsigned(c[0]) + c[1] + c[2] + c[3];
        if (sum >= cap) {
            for (auto& v : c) v = static_cast<std::uint16_t>((v + 1u) >> 1);
        }
    }
    void predict(Distribution& out, unsigned smooth) const {
        const std::uint64_t t = std::uint64_t(c[0]) + c[1] + c[2] + c[3];
        const std::uint64_t den = t * smooth + 4;
        const std::uint64_t rec = (UINT64_C(1) << 48) / den;
        for (unsigned b = 0; b < 4; ++b)
            out[b] = static_cast<std::uint32_t>(((std::uint64_t(c[b]) * smooth + 1) * rec) >> 32);
    }
};

struct BinaryCounter {
    std::array<std::uint16_t, 2> c{};
    void observe(unsigned bit, unsigned cap = 256) {
        if (c[bit] < 65000u) ++c[bit];
        if (unsigned(c[0]) + c[1] >= cap) {
            c[0] = static_cast<std::uint16_t>((c[0] + 1u) >> 1);
            c[1] = static_cast<std::uint16_t>((c[1] + 1u) >> 1);
        }
    }
    std::uint16_t predict(unsigned smooth = 24) const {
        const std::uint64_t t = std::uint64_t(c[0]) + c[1];
        return static_cast<std::uint16_t>(std::clamp<std::uint64_t>(
            ((std::uint64_t(c[1]) * smooth + 1) * 4096) / (t * smooth + 2), 1, 4095));
    }
};

class OrbitTableV3 {
    struct Row {
        std::uint64_t tag = UINT64_MAX;
        std::array<std::uint16_t, 4> counts{};
    };
    unsigned length_, ways_;
    std::uint64_t mask_;
    std::vector<Row> rows_;
    std::size_t index_mask_;

    Row* find(std::uint64_t key, bool insert) {
        const std::size_t first = (hash_word(key) & index_mask_) & ~std::size_t(ways_ - 1);
        unsigned found = ways_;
        for (unsigned i = 0; i < ways_; ++i) {
            if (rows_[first + i].tag == key) { found = i; break; }
        }
        if (found == ways_) {
            if (!insert) return nullptr;
            found = ways_ - 1;
            rows_[first + found] = {key, {}};
        }
        const Row row = rows_[first + found];
        for (unsigned i = found; i > 0; --i) rows_[first + i] = rows_[first + i - 1];
        rows_[first] = row;
        return &rows_[first];
    }

public:
    OrbitTableV3(unsigned length, unsigned ways, unsigned bits)
        : length_(length), ways_(ways), mask_(word_mask(length)),
          rows_(std::size_t{1} << bits), index_mask_(rows_.size() - 1) {}

    std::array<std::uint64_t, 3> counts(std::uint64_t word) {
        word &= mask_;
        const auto canon = canonical_word(word, length_);
        const Row* found = find(canon.word, false);
        if (!found) return {};
        const unsigned t = canon.transform;
        const auto r = reverse_word(word, length_);
        const bool rev_fix = (r == word);
        const bool rc_fix = ((r ^ mask_) == word);
        const std::uint64_t raw = found->counts[t];
        const std::uint64_t rc = rc_fix ? raw : raw + found->counts[t ^ 3];
        std::uint64_t pooled = 0;
        for (auto v : found->counts) pooled += v;
        if (rev_fix || rc_fix) pooled >>= 1;
        return {raw, rc, pooled};
    }

    void observe(std::uint64_t word) {
        word &= mask_;
        const auto canon = canonical_word(word, length_);
        Row& row = *find(canon.word, true);
        const unsigned t = canon.transform;
        if (row.counts[t] < 65000u) ++row.counts[t];
        const auto r = reverse_word(word, length_);
        if (r == word && row.counts[t ^ 1] < 65000u) ++row.counts[t ^ 1];
        else if ((r ^ mask_) == word && row.counts[t ^ 3] < 65000u) ++row.counts[t ^ 3];
        const unsigned sum = unsigned(row.counts[0]) + row.counts[1] + row.counts[2] + row.counts[3];
        if (sum >= 256u) {
            for (auto& v : row.counts) v = static_cast<std::uint16_t>((v + 1u) >> 1);
        }
    }
};

class MatchModelV3 {
    struct Slot {
        std::uint64_t tag = UINT64_MAX;
        std::array<std::uint32_t, 4> pos{};
    };
    struct Active {
        std::size_t next = 0, recent = 0;
        int dir = 1;
        unsigned comp = 0, len = 0, err = 0, base = 0;
        std::uint32_t win = 0;
    };
    std::array<std::vector<Slot>, 4> slots_;
    std::array<Active, 4> active_{};
    unsigned active_n_ = 0;

    static unsigned pop(std::uint32_t w) {
        unsigned n = 0;
        while (w) { w &= w - 1; ++n; }
        return n;
    }

public:
    MatchModelV3() {
        slots_[0].resize(std::size_t{1} << 21);
        slots_[1].resize(std::size_t{1} << 22);
        slots_[2].resize(std::size_t{1} << 22);
        slots_[3].resize(std::size_t{1} << 22);
    }
    void reset() { active_n_ = 0; }

    Distribution predict(const Bytes& hist, std::size_t pos, std::uint64_t ctx, std::size_t since_reset, unsigned& out_len, unsigned k0 = 11u) {
        const std::array<unsigned, 4> ks{k0, 13u, 16u, 24u};
        std::array<Active, 64> choices{};
        unsigned count = 0;
        auto add = [&](Active c) {
            if (c.next >= pos) return;
            c.base = hist[c.next] ^ c.comp;
            for (unsigned j = 0; j < count; ++j) {
                if (choices[j].next == c.next && choices[j].dir == c.dir && choices[j].comp == c.comp) {
                    if (c.len > choices[j].len) choices[j] = c;
                    return;
                }
            }
            if (count < choices.size()) choices[count++] = c;
        };
        bool keep_top = (active_n_ == 4);
        for (unsigned i = 0; i < active_n_; ++i) {
            add(active_[i]);
            if (int(active_[i].len) - 4 * int(active_[i].err) <= 96) keep_top = false;
        }
        for (unsigned w = 0; w < 4; ++w) {
            const unsigned k = ks[w];
            if (since_reset < k) continue;
            const auto word = ctx & word_mask(k);
            const auto canon = canonical_word(word, k);
            auto& row = slots_[w][hash_word(canon.word) & (slots_[w].size() - 1)];
            if (row.tag == canon.word) {
                for (auto packed : row.pos) {
                    if (keep_top) break;
                    if (!packed) continue;
                    const auto prev = static_cast<std::size_t>(packed >> 2);
                    const auto g = canon.transform ^ static_cast<unsigned>(packed & 3);
                    const auto rev = reverse_word(word, k);
                    const unsigned stab = (rev == word) ? 1u : (((rev ^ word_mask(k)) == word) ? 3u : 0u);
                    for (unsigned alt = 0; alt < (stab ? 2u : 1u); ++alt) {
                        const auto tr = g ^ (alt ? stab : 0u);
                        const bool rev_dir = (tr & 1u) != 0;
                        if (rev_dir && prev <= k) continue;
                        Active c;
                        c.next = rev_dir ? prev - k - 1 : prev;
                        c.dir = rev_dir ? -1 : 1;
                        c.comp = (tr & 2u) ? 3u : 0u;
                        c.recent = prev;
                        if (c.next >= pos) continue;
                        if (transform_word(row.tag, k, canon.transform) != word) continue;
                        c.len = k;
                        for (unsigned d = k + 1; d <= 512 && d <= since_reset; ++d) {
                            std::size_t src = 0;
                            if (rev_dir) { if (d >= pos - c.next) break; src = c.next + d; }
                            else { if (d > c.next) break; src = c.next - d; }
                            if ((hist[src] ^ c.comp) != hist[pos - d]) break;
                            ++c.len;
                        }
                        add(c);
                    }
                }
            } else {
                row.tag = canon.word;
                row.pos = {};
            }
            if (pos < (std::size_t{1} << 30)) {
                row.pos[3] = row.pos[2];
                row.pos[2] = row.pos[1];
                row.pos[1] = row.pos[0];
                row.pos[0] = static_cast<std::uint32_t>((pos << 2) | canon.transform);
            }
        }
        auto better = [](const Active& a, const Active& b) {
            const int sa = int(a.len) - 4 * int(a.err), sb = int(b.len) - 4 * int(b.err);
            if (sa != sb) return sa > sb;
            return a.recent > b.recent;
        };
        active_n_ = std::min(4u, count);
        for (unsigned i = 0; i < active_n_; ++i) {
            unsigned best = i;
            for (unsigned j = i + 1; j < count; ++j) if (better(choices[j], choices[best])) best = j;
            std::swap(choices[i], choices[best]);
            active_[i] = choices[i];
        }
        if (!active_n_) { out_len = 0; return {16384, 16384, 16384, 16384}; }
        out_len = active_[0].len;
        std::array<std::uint64_t, 4> votes{};
        std::uint64_t total = 0;
        for (unsigned i = 0; i < active_n_; ++i) {
            const auto& c = active_[i];
            const unsigned str = std::max(1u, c.len - std::min(c.len, 9u));
            const unsigned wt = std::max(1u, str * str / (1u + 6u * c.err));
            const unsigned err_p = std::max(4u, 4096u / (1u + str * 2u));
            for (unsigned b = 0; b < 4; ++b)
                votes[b] += std::uint64_t(wt) * (b == c.base ? scale - 3u * err_p : err_p);
            total += wt;
        }
        Distribution res{};
        for (unsigned b = 0; b < 4; ++b) res[b] = static_cast<std::uint32_t>(votes[b] / total);
        return res;
    }

    void observe(unsigned base) {
        unsigned kept = 0;
        for (unsigned i = 0; i < active_n_; ++i) {
            auto c = active_[i];
            const bool mis = (base != c.base);
            c.win = (c.win << 1) | unsigned(mis);
            c.err = pop(c.win);
            if (c.err > 2u) continue;
            if (c.dir < 0 && c.next == 0) continue;
            c.next = c.dir > 0 ? c.next + 1 : c.next - 1;
            c.len = std::min(4096u, c.len + 1);
            active_[kept++] = c;
        }
        active_n_ = kept;
    }
};

class PredictorV3 {
    // Experts (37 4-way): dense 0..9, orbits 10..24, codon 25..28, sparse 29..31, run 32, match 33, rc 34..36
    static constexpr unsigned NUM_EXPERTS = 37;
    static constexpr std::array<unsigned, 5> orbit_orders{10, 12, 15, 18, 22};

    std::vector<CompactCounter> d0_, d1_, d2_, d3_, d4_, d5_, d6_, d7_, d8_, d9_;
    std::vector<CompactCounter> codon2_, codon5_, codon8_, codon_wob_;
    std::vector<CompactCounter> skip1_, skip2_, gapped_, run_tbl_;
    std::vector<std::uint16_t> rc8_, rc9_, rc10_;
    std::vector<OrbitTableV3> orbits_;
    MatchModelV3 match_;
    using Weights = std::array<std::int32_t, NUM_EXPERTS>;
    Weights global_{};
    std::array<Weights, 96> local_{};
    std::array<Distribution, NUM_EXPERTS> preds_{};

    DnaStrategy strategy_ = DnaStrategy::coding_large;
    unsigned bucket_ = 0;
    unsigned last_apm_ = 0;
    std::vector<std::array<std::uint16_t, 4>> apm_;
    std::uint64_t ctx_ = 0;
    std::size_t since_reset_ = 0;
    unsigned last_base_ = 0, run_len_ = 0;

    void update_weights(Weights& w, unsigned base, int decay, int clamp_min) {
        const auto& lt = weight_tables();
        std::int32_t mx = std::numeric_limits<std::int32_t>::min();
        for (unsigned e = 0; e < NUM_EXPERTS; ++e) {
            const unsigned p = std::max(1u, std::min(65535u, preds_[e][base]));
            w[e] = (w[e] * decay) / 65536 + lt.log_probability[p];
            mx = std::max(mx, w[e]);
        }
        for (auto& v : w) v = std::max(clamp_min, v - mx);
    }
    static void eval_rc_exact(const std::vector<std::uint16_t>& tbl, std::uint64_t ctx, unsigned len, unsigned smooth, Distribution& out) {
        const std::uint64_t m = word_mask(len); std::array<std::uint64_t, 4> c{};
        for (unsigned b = 0; b < 4; ++b) { const auto w = ((ctx << 2) | b) & m; const auto rc = reverse_word(w, len) ^ m; c[b] = (w == rc) ? tbl[w] : (std::uint64_t(tbl[w]) + tbl[rc]); }
        const std::uint64_t rec = (UINT64_C(1) << 48) / ((c[0] + c[1] + c[2] + c[3]) * smooth + 4);
        for (unsigned b = 0; b < 4; ++b) out[b] = static_cast<std::uint32_t>(((c[b] * smooth + 1) * rec) >> 32);
    }
    static void obs_rc_exact(std::vector<std::uint16_t>& tbl, std::uint64_t word, unsigned len, unsigned cap) {
        const std::uint64_t w = word & word_mask(len); if (tbl[w] < 65000u) ++tbl[w];
        const std::uint64_t bc = w & ~std::uint64_t(3u);
        if (unsigned(tbl[bc]) + tbl[bc | 1] + tbl[bc | 2] + tbl[bc | 3] >= cap)
            for (unsigned b = 0; b < 4; ++b) tbl[bc | b] = static_cast<std::uint16_t>((tbl[bc | b] + 1u) >> 1);
    }

    static void eval_orbit(OrbitTableV3& tbl, std::uint64_t ctx, unsigned smooth, Distribution* out3) {
        std::array<std::array<std::uint64_t, 4>, 3> ev{};
        for (unsigned b = 0; b < 4; ++b) {
            const auto c = tbl.counts((ctx << 2) | b);
            for (unsigned k = 0; k < 3; ++k) ev[k][b] = c[k];
        }
        for (unsigned k = 0; k < 3; ++k) {
            const auto& c = ev[k];
            const std::uint64_t tot = c[0] + c[1] + c[2] + c[3];
            const std::uint64_t rec = (UINT64_C(1) << 48) / (tot * smooth + 4);
            for (unsigned b = 0; b < 4; ++b)
                out3[k][b] = static_cast<std::uint32_t>(((c[b] * smooth + 1) * rec) >> 32);
        }
    }

    unsigned wobble_idx() const {
        if (since_reset_ < 9) return 0u;
        const unsigned ph = static_cast<unsigned>(since_reset_ % 3u);
        const unsigned b0 = unsigned(ctx_ & 3u);
        const unsigned b1 = unsigned((ctx_ >> 2) & 3u);
        const unsigned b2 = unsigned((ctx_ >> 4) & 3u);
        const unsigned b5 = unsigned((ctx_ >> 10) & 3u);
        const unsigned b8 = unsigned((ctx_ >> 16) & 3u);
        return ph * 1024u + ((b0 << 8) | (b1 << 6) | (b2 << 4) | (b5 << 2) | b8);
    }

public:
    explicit PredictorV3(DnaStrategy strategy = DnaStrategy::coding_large)
        : d0_(1), d1_(4), d2_(16), d3_(64), d4_(256), d5_(1024), d6_(4096), d7_(16384), d8_(65536), d9_(262144),
          codon2_(3 * 16), codon5_(3 * 1024), codon8_(3 * 65536), codon_wob_(3 * 1024),
          skip1_(256), skip2_(256), gapped_(1024), run_tbl_(4 * 8), rc8_(1u << 16, 0), rc9_(1u << 18, 0), rc10_(1u << 20, 0),
          strategy_(strategy),
          apm_(strategy == DnaStrategy::coding_large ? 3072u : 768u, {1024, 1024, 1024, 1024}) {
        orbits_.reserve(orbit_orders.size());
        orbits_.emplace_back(11, 4, 21);
        orbits_.emplace_back(13, 4, 21);
        orbits_.emplace_back(16, 4, 21);
        orbits_.emplace_back(19, 4, 21);
        orbits_.emplace_back(23, 4, 21);
    }

    void reset() {
        ctx_ = 0; since_reset_ = 0; last_base_ = 0; run_len_ = 0;
        match_.reset();
    }

    Cdf predict(const Bytes& hist, std::size_t pos) {
        d0_[0].predict(preds_[0], 32);
        d1_[since_reset_ >= 1 ? (ctx_ & 3u) : 0].predict(preds_[1], 24);
        d2_[since_reset_ >= 2 ? (ctx_ & 15u) : 0].predict(preds_[2], 24);
        d3_[since_reset_ >= 3 ? (ctx_ & 63u) : 0].predict(preds_[3], 24);
        d4_[since_reset_ >= 4 ? (ctx_ & 255u) : 0].predict(preds_[4], 24);
        d5_[since_reset_ >= 5 ? (ctx_ & 1023u) : 0].predict(preds_[5], 24);
        d6_[since_reset_ >= 6 ? (ctx_ & 4095u) : 0].predict(preds_[6], 20);
        d7_[since_reset_ >= 7 ? (ctx_ & 16383u) : 0].predict(preds_[7], 16);
        d8_[since_reset_ >= 8 ? (ctx_ & 65535u) : 0].predict(preds_[8], 12);
        d9_[since_reset_ >= 9 ? (ctx_ & 262143u) : 0].predict(preds_[9], 8);

        unsigned best_ord = 0;
        for (unsigned i = 0; i < orbit_orders.size(); ++i) {
            if (since_reset_ >= orbit_orders[i]) {
                eval_orbit(orbits_[i], ctx_, 20u + i * 8u, &preds_[10 + i * 3]);
                if (preds_[10 + i * 3] != Distribution{16384, 16384, 16384, 16384}) best_ord = i + 1;
            } else {
                for (unsigned k = 0; k < 3; ++k) preds_[10 + i * 3 + k] = {16384, 16384, 16384, 16384};
            }
        }
        const unsigned ph = static_cast<unsigned>(since_reset_ % 3u);
        codon2_[ph * 16u + unsigned(ctx_ & 15u)].predict(preds_[25], 24);
        codon5_[ph * 1024u + unsigned(ctx_ & 1023u)].predict(preds_[26], 20);
        codon8_[ph * 65536u + (since_reset_ >= 8 ? unsigned(ctx_ & 65535u) : 0u)].predict(preds_[27], 12);
        codon_wob_[wobble_idx()].predict(preds_[28], 16);

        const unsigned s1 = since_reset_ >= 6 ? unsigned(((ctx_ & 3u) << 6) | (((ctx_ >> 4) & 3u) << 4) | (((ctx_ >> 6) & 3u) << 2) | ((ctx_ >> 10) & 3u)) : 0u;
        const unsigned s2 = since_reset_ >= 6 ? unsigned((((ctx_ >> 2) & 3u) << 6) | (((ctx_ >> 4) & 3u) << 4) | (((ctx_ >> 8) & 3u) << 2) | ((ctx_ >> 10) & 3u)) : 0u;
        const unsigned g3 = since_reset_ >= 7 ? unsigned(((ctx_ & 15u) << 6) | (((ctx_ >> 6) & 15u) << 2) | ((ctx_ >> 12) & 3u)) : 0u;
        skip1_[s1].predict(preds_[29], 24);
        skip2_[s2].predict(preds_[30], 24);
        gapped_[g3].predict(preds_[31], 20);
        run_tbl_[last_base_ * 8u + std::min(7u, run_len_)].predict(preds_[32], strategy_ == DnaStrategy::low_complex ? 48u : 24u);

        unsigned match_len = 0;
        const unsigned k0 = (strategy_ == DnaStrategy::short_reads) ? 9u : 11u;
        preds_[33] = match_.predict(hist, pos, ctx_, since_reset_, match_len, k0);
        if (since_reset_ >= 7) eval_rc_exact(rc8_, ctx_, 8, 16, preds_[34]); else preds_[34] = {16384, 16384, 16384, 16384};
        if (since_reset_ >= 8) eval_rc_exact(rc9_, ctx_, 9, 12, preds_[35]); else preds_[35] = {16384, 16384, 16384, 16384};
        if (since_reset_ >= 9) eval_rc_exact(rc10_, ctx_, 10, 8, preds_[36]); else preds_[36] = {16384, 16384, 16384, 16384};
        const unsigned match_bin = match_len == 0 ? 0u : (match_len < 24u ? 1u : 2u);
        bucket_ = ((best_ord * 3u + match_bin) * 3u + ph) % 96u;

        // 4-ary Log-linear + linear hybrid
        const auto& lt = weight_tables();
        std::array<std::int64_t, 4> log_mix{};
        std::array<std::uint64_t, 4> lin_mix{};
        std::uint64_t lsum = 0;
        for (unsigned e = 0; e < NUM_EXPERTS; ++e) {
            if (e > 0 && preds_[e][0] == 16384u && preds_[e][1] == 16384u && preds_[e][2] == 16384u) continue;
            const std::uint32_t wg = lt.weight[static_cast<std::size_t>(-global_[e])];
            const std::uint32_t wl = lt.weight[static_cast<std::size_t>(-local_[bucket_][e])];
            const std::uint64_t w = (std::uint64_t(wg) + std::uint64_t(wl) * 2u) >> 10;
            lsum += w;
            for (unsigned b = 0; b < 4; ++b) {
                const unsigned p = std::max(1u, std::min(65535u, preds_[e][b]));
                log_mix[b] += std::int64_t(w) * std::int64_t(lt.log_probability[p]);
                lin_mix[b] += w * std::uint64_t(p);
            }
        }
        std::array<std::int32_t, 4> avg_log{};
        std::int32_t max_log = std::numeric_limits<std::int32_t>::min();
        for (unsigned b = 0; b < 4; ++b) {
            avg_log[b] = static_cast<std::int32_t>(log_mix[b] / static_cast<std::int64_t>(std::max<std::uint64_t>(1, lsum)));
            max_log = std::max(max_log, avg_log[b]);
        }
        std::array<std::uint64_t, 4> mixed{};
        std::uint64_t sum = 0;
        unsigned top = 0;
        for (unsigned b = 0; b < 4; ++b) {
            const std::int32_t diff = std::min<std::int32_t>(8192, max_log - avg_log[b]);
            const std::uint64_t geo = lt.weight[static_cast<std::size_t>(diff)] >> 8;
            const std::uint64_t lin = lin_mix[b] / std::max<std::uint64_t>(1, lsum);
            mixed[b] = geo * 4u + lin;
            if (mixed[b] > mixed[top]) top = b;
            sum += mixed[b];
        }
        const unsigned top_p = static_cast<unsigned>(std::min<std::uint64_t>(15u, (mixed[top] * 16u) / std::max<std::uint64_t>(1, sum)));
        const unsigned apm_ctx_base = (strategy_ == DnaStrategy::coding_large && since_reset_ >= 2) ? unsigned(ctx_ & 15u) : last_base_;
        last_apm_ = (((apm_ctx_base * 3u + ph) * 4u + top) * 16u + top_p) % static_cast<unsigned>(apm_.size());
        sum = 0;
        for (unsigned b = 0; b < 4; ++b) {
            mixed[b] = std::max<std::uint64_t>(1u, (mixed[b] * apm_[last_apm_][b]) >> 10);
            sum += mixed[b];
        }
        Cdf cdf{0, 0, 0, 0, scale};
        std::uint64_t cum = 0;
        for (unsigned b = 0; b < 3; ++b) {
            cum += mixed[b];
            cdf[b + 1] = b + 1 + static_cast<std::uint32_t>(cum * (scale - 4) / std::max<std::uint64_t>(1, sum));
        }
        return cdf;
    }

    void observe(unsigned base, std::size_t) {
        const int ldecay = (strategy_ == DnaStrategy::short_reads || strategy_ == DnaStrategy::low_complex) ? 64200 : 64800;
        update_weights(global_, base, 65280, -8192);
        update_weights(local_[bucket_], base, ldecay, -8192);
        auto& ac = apm_[last_apm_];
        for (unsigned b = 0; b < 4; ++b) {
            const int target = (b == base) ? 1536 : 853;
            ac[b] = static_cast<std::uint16_t>(std::clamp(int(ac[b]) + ((target - int(ac[b])) >> 6), 512, 2048));
        }
        match_.observe(base);

        d0_[0].observe(base, 512);
        if (since_reset_ >= 1) d1_[ctx_ & 3u].observe(base, 256);
        if (since_reset_ >= 2) d2_[ctx_ & 15u].observe(base, 256);
        if (since_reset_ >= 3) d3_[ctx_ & 63u].observe(base, 256);
        if (since_reset_ >= 4) d4_[ctx_ & 255u].observe(base, 192);
        if (since_reset_ >= 5) d5_[ctx_ & 1023u].observe(base, 128);
        if (since_reset_ >= 6) d6_[ctx_ & 4095u].observe(base, 128);
        if (since_reset_ >= 7) d7_[ctx_ & 16383u].observe(base, 96);
        if (since_reset_ >= 8) d8_[ctx_ & 65535u].observe(base, 96);
        if (since_reset_ >= 9) d9_[ctx_ & 262143u].observe(base, 64);

        const unsigned ph = static_cast<unsigned>(since_reset_ % 3u);
        codon2_[ph * 16u + unsigned(ctx_ & 15u)].observe(base, 192);
        codon5_[ph * 1024u + unsigned(ctx_ & 1023u)].observe(base, 128);
        if (since_reset_ >= 8) codon8_[ph * 65536u + unsigned(ctx_ & 65535u)].observe(base, 96);
        if (since_reset_ >= 9) codon_wob_[wobble_idx()].observe(base, 128);

        if (since_reset_ >= 6) {
            const unsigned s1 = unsigned(((ctx_ & 3u) << 6) | (((ctx_ >> 4) & 3u) << 4) | (((ctx_ >> 6) & 3u) << 2) | ((ctx_ >> 10) & 3u));
            const unsigned s2 = unsigned((((ctx_ >> 2) & 3u) << 6) | (((ctx_ >> 4) & 3u) << 4) | (((ctx_ >> 8) & 3u) << 2) | ((ctx_ >> 10) & 3u));
            skip1_[s1].observe(base, 128);
            skip2_[s2].observe(base, 128);
        }
        if (since_reset_ >= 7) {
            const unsigned g3 = unsigned(((ctx_ & 15u) << 6) | (((ctx_ >> 6) & 15u) << 2) | ((ctx_ >> 12) & 3u));
            gapped_[g3].observe(base, 128);
        }
        run_tbl_[last_base_ * 8u + std::min(7u, run_len_)].observe(base, 128);
        if (since_reset_ > 0 && base == last_base_) ++run_len_;
        else { last_base_ = base; run_len_ = 1; }

        const std::uint64_t word = (ctx_ << 2) | base;
        if (since_reset_ >= 7) obs_rc_exact(rc8_, word, 8, 96);
        if (since_reset_ >= 8) obs_rc_exact(rc9_, word, 9, 96);
        if (since_reset_ >= 9) obs_rc_exact(rc10_, word, 10, 64);
        for (unsigned i = 0; i < orbit_orders.size(); ++i)
            if (since_reset_ >= orbit_orders[i]) orbits_[i].observe(word);
        ctx_ = word & word_mask(24);
        ++since_reset_;
    }
};

} // namespace v4zip::detail
