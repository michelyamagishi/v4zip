#pragma once
#include "range.hpp"
#include "weights.hpp"
#include <algorithm>
#include <array>
#include <cstdint>
#include <vector>

namespace v4zip::detail {

// Multi-expert adaptive range coder for FASTQ quality strings
class QualityModel {
    static constexpr unsigned SYMS = 96; // ASCII 33..126 + fallback
    static constexpr unsigned NUM_EXP = 7;

    struct HashRow {
        std::uint32_t tag = UINT32_MAX;
        std::uint16_t total = 0;
        std::array<std::uint16_t, SYMS> c{};
    };
    struct DenseRow {
        std::uint16_t total = 0;
        std::array<std::uint16_t, SYMS> c{};
    };

    DenseRow o0_{};
    std::vector<DenseRow> o1_;
    std::vector<DenseRow> o2_direct_; // 64 x 64 exact for p1, p2 < 64
    std::vector<HashRow> o2_hash_;    // fallback for >= 64
    std::vector<HashRow> o3_;         // 2-way associative (16384 rows = 8192 buckets)
    std::vector<DenseRow> pos_p1_;    // 128 x 64 exact pos x p1
    std::vector<DenseRow> run_tbl_;   // 8 x 64
    std::vector<DenseRow> mean_tbl_;  // 16 x 32

    std::array<std::int32_t, NUM_EXP> gw_{};
    std::array<std::array<std::int32_t, NUM_EXP>, 16> lw_{};
    std::array<std::array<std::uint32_t, SYMS>, NUM_EXP> ep_{};
    unsigned last_pb_ = 0;

    static unsigned sym_of(std::uint8_t c) {
        return (c >= 33 && c <= 126) ? unsigned(c - 33) : 95u;
    }

    static std::uint32_t mix_hash(std::uint32_t x) {
        x ^= x >> 16; x *= 0x7feb352du;
        x ^= x >> 15; x *= 0x846ca68bu;
        return x ^ (x >> 16);
    }

    static void row_to_probs(const std::array<std::uint16_t, SYMS>& c, std::uint16_t total, unsigned smooth, std::array<std::uint32_t, SYMS>& out) {
        const std::uint64_t den = std::uint64_t(total) * smooth + SYMS;
        const std::uint64_t rec = (UINT64_C(1) << 44) / den;
        for (unsigned s = 0; s < SYMS; ++s) {
            out[s] = static_cast<std::uint32_t>(std::clamp<std::uint64_t>(((std::uint64_t(c[s]) * smooth + 1u) * rec) >> 28, 1u, 65535u));
        }
    }

    static void bump_dense(DenseRow& row, unsigned sym, std::uint16_t cap) {
        if (row.c[sym] < 60000u) { ++row.c[sym]; ++row.total; }
        if (row.total >= cap) {
            std::uint16_t t = 0;
            for (auto& v : row.c) { v = static_cast<std::uint16_t>(v >> 1); t = static_cast<std::uint16_t>(t + v); }
            row.total = t;
        }
    }

    static const HashRow* find_2way(const std::vector<HashRow>& tbl, std::uint32_t tag) {
        const std::size_t base = (mix_hash(tag) & ((tbl.size() >> 1) - 1)) << 1;
        if (tbl[base].tag == tag) return &tbl[base];
        if (tbl[base + 1].tag == tag) return &tbl[base + 1];
        return nullptr;
    }

    static void bump_2way(std::vector<HashRow>& tbl, std::uint32_t tag, unsigned sym, std::uint16_t cap) {
        const std::size_t base = (mix_hash(tag) & ((tbl.size() >> 1) - 1)) << 1;
        std::size_t slot = base;
        if (tbl[base].tag == tag) slot = base;
        else if (tbl[base + 1].tag == tag) slot = base + 1;
        else {
            slot = (tbl[base].total <= tbl[base + 1].total) ? base : (base + 1);
            tbl[slot].tag = tag;
            tbl[slot].total = 0;
            tbl[slot].c.fill(0);
        }
        auto& row = tbl[slot];
        if (row.c[sym] < 60000u) { ++row.c[sym]; ++row.total; }
        if (row.total >= cap) {
            std::uint16_t t = 0;
            for (auto& v : row.c) { v = static_cast<std::uint16_t>(v >> 1); t = static_cast<std::uint16_t>(t + v); }
            row.total = t;
        }
    }

public:
    QualityModel()
        : o1_(SYMS), o2_direct_(64 * 64), o2_hash_(512), o3_(16384),
          pos_p1_(128 * 64), run_tbl_(8 * 64), mean_tbl_(16 * 32) {}

    void build_cdf(unsigned p1, unsigned p2, unsigned p3, unsigned pos, unsigned run_len, unsigned mean_q,
                   std::array<std::uint32_t, SYMS + 1>& cdf) {
        const std::uint32_t t2 = (p1 << 8) | p2;
        const std::uint32_t t3 = (p1 << 16) | (p2 << 8) | p3;
        const unsigned pp_idx = std::min(127u, pos) * 64u + std::min(63u, p1);
        const unsigned run_idx = std::min(7u, run_len) * 64u + std::min(63u, p1);
        const unsigned mean_idx = std::min(15u, mean_q / 3u) * 32u + std::min(31u, pos >> 1);

        row_to_probs(o0_.c, o0_.total, 16, ep_[0]);
        row_to_probs(o1_[p1].c, o1_[p1].total, 24, ep_[1]);
        if (p1 < 64u && p2 < 64u) {
            const auto& r2 = o2_direct_[p1 * 64u + p2];
            row_to_probs(r2.c, r2.total, 32, ep_[2]);
        } else if (const auto* r2 = find_2way(o2_hash_, t2)) {
            row_to_probs(r2->c, r2->total, 32, ep_[2]);
        } else {
            ep_[2].fill(65536u / SYMS);
        }
        if (const auto* r3 = find_2way(o3_, t3)) {
            row_to_probs(r3->c, r3->total, 40, ep_[3]);
        } else {
            ep_[3].fill(65536u / SYMS);
        }
        row_to_probs(pos_p1_[pp_idx].c, pos_p1_[pp_idx].total, 28, ep_[4]);
        row_to_probs(run_tbl_[run_idx].c, run_tbl_[run_idx].total, 24, ep_[5]);
        row_to_probs(mean_tbl_[mean_idx].c, mean_tbl_[mean_idx].total, 24, ep_[6]);

        last_pb_ = std::min(15u, pos >> 2);
        const auto& lt = weight_tables();
        std::array<std::int64_t, SYMS> log_mix{};
        std::array<std::uint64_t, SYMS> lin_mix{};
        std::uint64_t lsum = 0;
        for (unsigned e = 0; e < NUM_EXP; ++e) {
            const std::uint32_t wg = lt.weight[static_cast<std::size_t>(std::clamp(-gw_[e], 0, 8192))];
            const std::uint32_t wl = lt.weight[static_cast<std::size_t>(std::clamp(-lw_[last_pb_][e], 0, 8192))];
            const std::uint64_t w = (std::uint64_t(wg) + std::uint64_t(wl) * 2u) >> 10;
            lsum += w;
            for (unsigned s = 0; s < SYMS; ++s) {
                const unsigned p = ep_[e][s];
                log_mix[s] += std::int64_t(w) * std::int64_t(lt.log_probability[p]);
                lin_mix[s] += w * std::uint64_t(p);
            }
        }
        std::array<std::int32_t, SYMS> avg_log{};
        std::int32_t max_log = std::numeric_limits<std::int32_t>::min();
        for (unsigned s = 0; s < SYMS; ++s) {
            avg_log[s] = static_cast<std::int32_t>(log_mix[s] / static_cast<std::int64_t>(std::max<std::uint64_t>(1, lsum)));
            max_log = std::max(max_log, avg_log[s]);
        }
        std::array<std::uint64_t, SYMS> mixed{};
        std::uint64_t sum = 0;
        for (unsigned s = 0; s < SYMS; ++s) {
            const std::int32_t diff = std::min<std::int32_t>(8192, max_log - avg_log[s]);
            const std::uint64_t geo = lt.weight[static_cast<std::size_t>(diff)] >> 8;
            const std::uint64_t lin = lin_mix[s] / std::max<std::uint64_t>(1, lsum);
            mixed[s] = geo * 3u + lin;
            sum += mixed[s];
        }
        cdf[0] = 0;
        cdf[SYMS] = scale;
        std::uint64_t cum = 0;
        for (unsigned s = 0; s + 1 < SYMS; ++s) {
            cum += mixed[s];
            cdf[s + 1] = s + 1 + static_cast<std::uint32_t>(cum * (scale - SYMS) / std::max<std::uint64_t>(1, sum));
        }
    }

    void observe(unsigned p1, unsigned p2, unsigned p3, unsigned pos, unsigned run_len, unsigned mean_q, unsigned s) {
        const auto& lt = weight_tables();
        auto upd_w = [&](std::array<std::int32_t, NUM_EXP>& w, int decay) {
            std::int32_t mx = std::numeric_limits<std::int32_t>::min();
            for (unsigned e = 0; e < NUM_EXP; ++e) {
                w[e] = (w[e] * decay) / 65536 + lt.log_probability[ep_[e][s]];
                mx = std::max(mx, w[e]);
            }
            for (auto& v : w) v = std::max(-8192, v - mx);
        };
        upd_w(gw_, 65280);
        upd_w(lw_[last_pb_], 64800);

        const std::uint32_t t2 = (p1 << 8) | p2;
        const std::uint32_t t3 = (p1 << 16) | (p2 << 8) | p3;
        const unsigned pp_idx = std::min(127u, pos) * 64u + std::min(63u, p1);
        const unsigned run_idx = std::min(7u, run_len) * 64u + std::min(63u, p1);
        const unsigned mean_idx = std::min(15u, mean_q / 3u) * 32u + std::min(31u, pos >> 1);

        bump_dense(o0_, s, 2048u);
        bump_dense(o1_[p1], s, 512u);
        if (p1 < 64u && p2 < 64u) bump_dense(o2_direct_[p1 * 64u + p2], s, 512u);
        else bump_2way(o2_hash_, t2, s, 384u);
        bump_2way(o3_, t3, s, 384u);
        bump_dense(pos_p1_[pp_idx], s, 512u);
        bump_dense(run_tbl_[run_idx], s, 256u);
        bump_dense(mean_tbl_[mean_idx], s, 256u);
    }

    static Bytes encode(const Bytes& qual, const std::vector<std::size_t>& read_starts, std::size_t) {
        for (auto c : qual) if (c < 33 || c > 126) return {};
        QualityModel m;
        RangeEncoder enc;
        std::array<std::uint32_t, SYMS + 1> cdf{};
        std::size_t b = 0, read_pos = 0, q_sum = 0;
        unsigned p1 = 40, p2 = 40, p3 = 40, run_len = 0;
        for (std::size_t i = 0; i < qual.size(); ++i) {
            while (b < read_starts.size() && read_starts[b] == i) {
                p1 = 40; p2 = 40; p3 = 40; read_pos = 0; run_len = 0; q_sum = 0; ++b;
            }
            const unsigned mean_q = read_pos ? static_cast<unsigned>(q_sum / read_pos) : 40u;
            const unsigned s = sym_of(qual[i]);
            m.build_cdf(p1, p2, p3, static_cast<unsigned>(read_pos), run_len, mean_q, cdf);
            enc.put(s, cdf);
            m.observe(p1, p2, p3, static_cast<unsigned>(read_pos), run_len, mean_q, s);
            run_len = (read_pos > 0 && s == p1) ? (run_len + 1) : 1u;
            p3 = p2; p2 = p1; p1 = s; q_sum += s; ++read_pos;
        }
        return enc.finish();
    }

    static Bytes decode(const Bytes& encoded, std::size_t expected, const std::vector<std::size_t>& read_starts) {
        QualityModel m;
        RangeDecoder dec(encoded);
        Bytes out;
        out.reserve(expected);
        std::array<std::uint32_t, SYMS + 1> cdf{};
        std::size_t b = 0, read_pos = 0, q_sum = 0;
        unsigned p1 = 40, p2 = 40, p3 = 40, run_len = 0;
        for (std::size_t i = 0; i < expected; ++i) {
            while (b < read_starts.size() && read_starts[b] == i) {
                p1 = 40; p2 = 40; p3 = 40; read_pos = 0; run_len = 0; q_sum = 0; ++b;
            }
            const unsigned mean_q = read_pos ? static_cast<unsigned>(q_sum / read_pos) : 40u;
            m.build_cdf(p1, p2, p3, static_cast<unsigned>(read_pos), run_len, mean_q, cdf);
            const unsigned s = dec.get(cdf);
            if (s >= 94u) fail("invalid quality symbol");
            out.push_back(static_cast<std::uint8_t>(s + 33));
            m.observe(p1, p2, p3, static_cast<unsigned>(read_pos), run_len, mean_q, s);
            run_len = (read_pos > 0 && s == p1) ? (run_len + 1) : 1u;
            p3 = p2; p2 = p1; p1 = s; q_sum += s; ++read_pos;
        }
        dec.finish();
        return out;
    }
};

} // namespace v4zip::detail
