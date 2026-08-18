# C++ Mastery → Quant Trading: Learning Roadmap

A personal roadmap covering modern C++ (memory management, moves, RAII),
concurrency & low-latency systems, and project ideas geared toward quant trading.

---

## Part 1 — Modern C++ fundamentals (moves, RAII, NRVO, smart pointers)

Going from `new`/`delete` to modern C++. The one reframing that ties it all
together: **you almost never write `new`/`delete` directly** — RAII objects
(containers, smart pointers) own memory and free it automatically, and *move
semantics* is just transferring that ownership cheaply instead of copying.

### Start here (free)
- **learncpp.com** — best free structured tutorial. Read the chapters on
  references / lvalues & rvalues, then "Move semantics and smart pointers".
- **cppreference.com** — authoritative reference (not a tutorial). Look up:
  - "Value categories" (lvalue / xvalue / prvalue) — foundation for moves
  - "Copy elision" (where RVO/NRVO are defined)
  - `std::move`, `std::forward`, `std::unique_ptr`, `std::shared_ptr`

### The book to get
- **"Effective Modern C++" — Scott Meyers.** *The* book for move semantics,
  perfect forwarding, RVO/NRVO, smart pointers, `auto`. ~40 short "Items".
- Companion overview: **"A Tour of C++" — Bjarne Stroustrup.**
- Deep dive: **"C++ Move Semantics: The Complete Guide" — Nicolai Josuttis.**

### Talks (YouTube — CppCon "Back to Basics" series)
- "Back to Basics: Move Semantics" (Klaus Iglberger) — clearest explanation
- "Back to Basics: RAII and the Rule of Zero/Five"
- "Back to Basics: Smart Pointers"
- "Back to Basics: Copy Elision / RVO"

### Concepts to internalize
1. **RAII** — resources owned by objects, freed in destructors. The core paradigm.
2. **Rule of Zero / Rule of Five** — when to write copy/move ctors & destructors.
3. **Smart pointers** — `unique_ptr` (single owner), `shared_ptr` (shared).
4. **Value categories** (lvalue vs rvalue) — explains *why* `std::move`/`&&` work.

### Quick `std::move` / return-by-value rules (from our discussion)
- `std::move` only helps for types that **own heap resources** (`string`,
  `vector`, `function`, `unique_ptr`, `Graph`). Pointless for `int`, `double`,
  `node`, pointers.
- **Constructor "sink" params:** take heavy args *by value*, then `std::move`
  into the member.
- **`return localVar;`** — plain return; NRVO/implicit move handle it.
  Never `return std::move(local)` (disables NRVO) and never return `T&&` of a
  local (dangling reference / UB).
- Don't `std::move` something you still use afterward (moved-from = valid but
  unspecified). `const` objects can't be moved (silently copy).
- Reference members (`const T&`) are **non-owning** — only safe if the referent
  outlives the object (fine for a long-lived `graph`, dangerous for a
  `std::function` typically passed as a temporary → own it by value instead).

---

## Part 2 — Concurrency, lock-free & low-latency systems

Leap from classic DSA: correctness is no longer just about the algorithm, it's
about the **memory model and ordering** — two threads can observe operations in
different orders unless you enforce otherwise.

### Concurrency & lock-free
- **"C++ Concurrency in Action" — Anthony Williams (2nd ed.)** ⭐ THE book.
  Threads → memory model → `std::atomic` → designing lock-free data structures
  → hazard pointers → memory reclamation.
- **"The Art of Multiprocessor Programming" — Herlihy & Shavit** ⭐ Definitive
  conceptual book: linearizability, lock-free/wait-free, ABA, concurrent
  structures. Java examples, language-agnostic ideas. Best fit given a strong
  DSA background.
- **"Concurrency with Modern C++" — Rainer Grimm** — practical on-ramp / newer
  stdlib features (`jthread`, coroutines, `atomic_ref`).

### Low-latency & performance
- **"Optimized C++" — Kurt Guntheroth** — practical perf tuning.
- **Agner Fog's optimization manuals (FREE, agner.org)** ⭐ microarchitecture,
  instruction timings — what people actually read for latency work.
- **"Systems Performance" — Brendan Gregg** — understand the machine (CPU,
  memory, scheduling, `perf`).

### Free papers (short, high-impact)
- **"What Every Programmer Should Know About Memory" — Ulrich Drepper** ⭐
  caches, memory hierarchy, false sharing. Read early.
- **"Is Parallel Programming Hard…" — Paul McKenney** (perfbook) — memory
  ordering, RCU, real-world lock-free. Advanced.

### Talks (search titles on YouTube)
- **Carl Cook — "When a Microsecond Is an Eternity"** (CppCon) ⭐ low-latency trading
- **Fedor Pikus — "Live Lock-Free or Deadlock"**, and his C++ atomics talks
- **David Gross — "Trading at Light Speed"**
- Timur Doumler / Anthony Williams — memory model talks

### Suggested order
1. Drepper's memory paper (background: caches, false sharing)
2. *C++ Concurrency in Action* Ch. 1–5 (threads, memory model, atomics)
3. *Art of Multiprocessor Programming* (linearizability + a few structures)
4. *C++ Concurrency in Action* Ch. 7 (build a lock-free stack/queue yourself)
5. Carl Cook's talk + Agner Fog (the low-latency mindset)

---

## Part 3 — Projects (learn modern C++ + build a quant portfolio)

### Tier 1 — Cement fundamentals (RAII, moves, templates)
1. **`Matrix` / linear-algebra class** — Rule of Five, moves, RVO, operator
   overloading, cache-friendly layout. (Later: expression templates.)
2. **Reimplement `unique_ptr` and `vector`** — the exercise that makes RAII and
   moves click permanently.
3. **`Decimal` / fixed-point type** — value semantics; prices aren't floats.

### Tier 2 — Quant-flavored building blocks
4. **Limit order book (LOB)** ⭐ the classic quant project. Price-time priority,
   add/cancel/execute. Then *optimize its latency* and tell that story.
5. **Market-data feed parser** — parse a binary (ITCH/FIX-style) protocol fast:
   zero-copy, `string_view`, no allocations in hot path.
6. **Options pricer: Black-Scholes + Monte Carlo** — numerics, `<random>`,
   template/`std::function` payoffs, parallelize with `thread`/`async`.

### Tier 3 — Systems & performance (what stands out)
7. **Lock-free SPSC ring buffer/queue** ⭐ `std::atomic`, memory ordering,
   false sharing, cache lines.
8. **Memory pool / arena allocator** — placement `new`, alignment, avoid heap
   in hot paths.
9. **Backtesting engine** — `Strategy` interface (pluggable-solver pattern),
   event loop, avoid lookahead bias, P&L / Sharpe.
10. **Latency-benchmarking harness** — nanosecond timing, p50/p99 histograms;
    tail latency matters more than average.

### Sequence
1. #2 (reimplement `vector`/`unique_ptr`) + #1 (`Matrix`) → moves/RAII lock in
2. #4 (order book) → signature quant project
3. #6 (options pricer) → numerics + parallelism
4. #7 (lock-free queue) + #8 (memory pool) → low-latency skills
5. #9 (backtester) → strategy + engineering together

### Portfolio force-multipliers
- **Always benchmark** (Google Benchmark). "Reduced order-book add latency
  800ns → 90ns by X" beats any feature.
- **Profile** (`perf`, `valgrind --tool=cachegrind`). Understand cache misses.
- **Test** (GoogleTest). Correctness under concurrency is where bugs hide.

> Single highest-leverage starter: **the limit order book (#4), then optimize
> and benchmark it.** Most direct "I can do quant C++" signal; pulls in data
> structures, moves, and latency thinking at once.
