# PlanarSeparator.cpp — audyt kroków 1–8, przypadki brzegowe i propozycja refaktoru

Dokument analizuje `cpp/separator/PlanarSeparator.cpp` (kroki 1–8 algorytmu
Liptona–Tarjana), wskazuje błędy poprawnościowe i przypadki brzegowe dla
zdegenerowanych grafów oraz proponuje lekki refaktor poprawiający czytelność.
Numeracja kroków jak w pracy i w kodzie.

Kontekst API (`PlanarGraphTools`):

- `planar_embedding_t = std::unordered_map<node, std::vector<node>>` — dostęp
  `embedding[v]` przez `operator[]` wstawia pustą listę dla wierzchołka bez
  krawędzi (to akurat wygodne, ale trzeba o tym pamiętać).
- `findPlanarEmbedding` **rzuca** `std::runtime_error("Graph have to be planar!")`
  dla grafu nieplanarnego.
- `makeMaximalPlanar` zwraca nowy graf (te same wierzchołki + dodane krawędzie
  triangulacji), wewnętrznie woła `make_biconnected_planar` + `make_maximal_planar`.

---

## 1. Błędy poprawnościowe (priorytet)

### 1.1. Wyszukiwanie `l0` — przedwczesny `break` (błąd rzędu separatora)

```cpp
NetworKit::node l0 = 0;
for (int i = 0; i <= l1; i++) {
  if (verticesAtLevel[i] + 2 * (l1 - i) <= 2 * sqrt(k)) {
    l0 = i;
  } else
    break;                     // <-- BŁĄD
}
```

Funkcja `f(i) = L(i) + 2*(l1 - i)` **nie jest monotoniczna** (składnik `2*(l1-i)`
maleje, ale `L(i)` skacze), więc zbiór poziomów spełniających warunek nie musi
być prefiksem. Konsekwencje:

- Jeśli `i = 0` nie spełnia warunku (np. `1 + 2*l1 > 2*sqrt(k)`), pętla przerywa
  od razu i zostaje `l0 = 0` z inicjalizacji — czyli **poziom, który nie spełnia
  warunku**. Wtedy gwarancja `|C| = O(√n)` przestaje obowiązywać.
- Jeśli spełniają `i = 0` i `i = 2`, a `i = 1` nie, kod zwróci `0` zamiast
  wyższego `2`.

Praca gwarantuje istnienie poprawnego `l0 ∈ [0, l1]` (dowód Th. 4), ale trzeba go
faktycznie znaleźć. Poprawka — szukaj od góry i weź pierwszy pasujący (najwyższy):

```cpp
NetworKit::node l0 = 0;
bool found0 = false;
for (int i = (int)l1; i >= 0; --i) {
  if (verticesAtLevel[i] + 2 * (l1 - i) <= 2 * std::sqrt((double)k)) {
    l0 = i; found0 = true; break;
  }
}
// found0 zawsze true na mocy dowodu; asercja pomaga wykryć regresje
assert(found0);
```

### 1.2. Wyszukiwanie `l2` — brak poziomu `r+1` i niepoprawny default

```cpp
NetworKit::node l2 = l1 + 1;
for (int i = l1 + 1; i < verticesAtLevel.size(); i++) {
  if (verticesAtLevel[i] + 2 * (i - l1 - 1) <= 2 * sqrt(n - k)) { l2 = i; break; }
}
```

- Wybór najmniejszego pasującego jest OK (praca chce „lowest level l2").
- Ale gdy **żaden** poziom w `[l1+1, maxLevel]` nie spełnia warunku, zostaje
  `l2 = l1 + 1` (default), który może nie spełniać warunku. Praca wprowadza
  dodatkowy **pusty poziom `r+1`**, który zawsze można wziąć (`L=0`). Pętla
  powinna rozważać też `i = maxLevel + 1` (z `L=0`), a nie zostawiać
  niepoprawnego defaultu.

Poprawka: dodaj wirtualny poziom `r+1` z `L = 0`:

```cpp
int maxLevel = (int)verticesAtLevel.size() - 1;
NetworKit::node l2 = maxLevel + 1;   // pusty poziom r+1 jako bezpieczny fallback
for (int i = l1 + 1; i <= maxLevel; i++) {
  if (verticesAtLevel[i] + 2 * (i - l1 - 1) <= 2 * std::sqrt((double)(n - k))) { l2 = i; break; }
}
// dla i = maxLevel+1: L=0, warunek 2*(maxLevel - l1) <= 2*sqrt(n-k) — sprawdź jawnie,
// a jeśli i to nie spełnia, i tak jest to poprawny separator poziomowy (patrz Lemat 3).
```

### 1.3. Graf pusty (n = 0) — crash

Dla `n = 0`:

- `areConnectedComponentsEligibleForPartition`: `largest = 0 < 2*0/3 = 0` → **false**,
  więc wchodzimy w gałąź główną (Step 3).
- `G` pusty, `root` pozostaje `NetworKit::none`, a
  `performBFSAndFindSpanningTree(G, none)` robi `vis[none] = true` → indeks poza
  zakresem → UB/crash.

Poprawka — wczesny guard na początku `run()`:

```cpp
if (graph.numberOfNodes() == 0) { hasRun = true; return; }
```

### 1.4. Wierzchołki spoza największej komponenty gubione w gałęzi głównej

W gałęzi `else` (największa komponenta ≥ 2n/3) pracujemy tylko na `G`
(największa komponenta) i wołamy `extractSeparatorAndPartitions(G, ...)`.
Wierzchołki **pozostałych komponent** nie trafiają do `A`, `B` ani `separator`.
Praca (Th. 4) wymaga rozszerzenia podziału na cały graf (pozostałe komponenty
rozdzielamy tak, by wyrównać koszty). To musi obsłużyć Step 10
(`extractSeparatorAndPartitions`) — do zaznaczenia jako wymóg, inaczej podział
nie pokrywa `V(G)`.

### 1.5. Normalizacja kosztu dla grafów niespójnych

W gałęzi głównej `n = G.numberOfNodes()` to rozmiar **największej komponenty**,
a progi (`n/2`, `sqrt(k)`, `sqrt(n-k)`) liczone są względem niej. Formalnie w
pracy, gdy komponenta ma koszt > 2/3 **całości**, poziom `l1` wyznaczamy tam,
gdzie skumulowany koszt osiąga 1/2 **całego grafu**, a nie połowy komponenty.
Dla grafu spójnego (`n = n_total`) jest OK; dla niespójnego z dominującą
komponentą progi względem rozmiaru komponenty mogą lekko rozjechać bilans.
Rekomendacja: progi liczyć względem całkowitej liczby wierzchołków grafu i/lub
udokumentować założenie.

### 1.6. Niezadeklarowana zmienna `side` (plik się nie kompiluje)

```cpp
extractSeparatorAndPartitions(G, lvl, l0, l2, cycle, x, side);  // 'side' nie istnieje
```

`side` powstaje dopiero w Step 9/10 (patrz `step9_ON_implementation.md`).
Do czasu implementacji Step 9 albo zaślepić (`std::vector<Side> side = ...`),
albo tymczasowo policzyć jednorazowym `markInsideOutside`.

---

## 2. Przypadki brzegowe / zdegenerowane grafy

| Wejście                                  | Zachowanie obecne                                                                                                    | Zalecenie                                                   |
| ---------------------------------------- | -------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------- |
| Graf pusty (n=0)                         | Crash (1.3)                                                                                                          | Wczesny `return`                                            |
| Graf nieplanarny                         | `findPlanarEmbedding` rzuca wyjątek (niełapany)                                                                      | Udokumentować/otoczyć try lub walidować                     |
| 1 wierzchołek                            | Gałąź główna → brak krawędzi nietrzewiowej → `fallbackLevelSeparator`                                                | OK, dodać test                                              |
| Sam las/drzewo                           | H po kontrakcji może nie mieć krawędzi nietrzewiowej do triangulacji; po `makeMaximalPlanar` OK; przy H<3 → fallback | Guard przed `makeMaximalPlanar` gdy `H.numberOfNodes() < 3` |
| Wszystkie wierzchołki izolowane          | `findSeparatorFromComponents` rozdziela                                                                              | OK, dodać test                                              |
| Niespójny, duża komponenta ≥ 2n/3        | Reszta komponent gubiona (1.4)                                                                                       | Rozszerzyć w Step 10                                        |
| `l0+1 == l2` (pusty środek)              | H = sam `x`, brak cyklu → fallback                                                                                   | OK                                                          |
| `n - k == 0` (k = n)                     | `sqrt(0)=0`; pętla `l2` pusta; default problematyczny (1.2)                                                          | Poziom `r+1`                                                |
| Self-loops / krawędzie równoległe        | Boost planarity może się wykrzaczyć                                                                                  | Założyć/wymusić graf prosty                                 |
| `makeMaximalPlanar` na 1–2 wierzchołkach | Ryzyko po stronie Boost                                                                                              | Guard j.w.                                                  |

Dodatkowe guardy sugerowane przed `makeMaximalPlanar`:

```cpp
if (H.numberOfNodes() < 3 || H.numberOfEdges() < 3) {
  fallbackLevelSeparator(lvl, l1);
  hasRun = true;
  return;
}
```

---

## 3. Drobne błędy / zapachy (bez zmiany logiki)

1. **`findNumberOfVerticesAtLevel`** deklaruje zwrot `std::vector<node>`, a zwraca
   licznik poziomów. `node`/`count` to ten sam typ bazowy, więc się kompiluje, ale
   sygnatura myli. Zmień typ zwracany na `std::vector<NetworKit::count>` i usuń
   komentarz „return type is wrong".
2. **Kopiowanie wektorów embeddingu w pętlach**: `auto neighborsV = embeddingH[v];`
   oraz `auto neighborsOfCycle1 = embeddingH[cycle[1]];` kopiują O(deg). Użyj
   `const auto& neighborsV = embeddingH.at(v);`.
3. **Trzykrotny `extractLargestConnectedComponent`** (w eligibility, w
   `findSeparatorFromComponents`, w gałęzi głównej) — kosztowne kopie tego samego.
   Policz raz i przekazuj.
4. **Brak `#include <numeric>`** dla `std::partial_sum` (działa tylko tranzytywnie).
   Dodaj jawnie.
5. **Porównania signed/unsigned**: `for (int i = ...; i < verticesAtLevel.size(); i++)`
   oraz `int i <= l1` (l1 to `node`/uint64). Ujednolić typy pętli.
6. **`std::find(...) == end()`** w detekcji orientacji cyklu nie jest zabezpieczone
   (zwróciłoby pozycję `size()`), choć teoretycznie sąsiedzi cyklu zawsze istnieją.
   Dodaj `assert` dla czytelności intencji.
7. **`inside`/`outside` jako `int`** — dla bardzo dużych grafów lepszy
   `NetworKit::count`/`int64_t` (koszt do n).

---

## 4. Propozycja lekkiego refaktoru (czytelność)

`run()` ma ~290 linii i miesza 8 kroków. Bez zmiany algorytmu warto rozbić gałąź
główną na prywatne metody/prosty kontekst. Proponowany podział:

```cpp
struct LevelStructure {
  std::vector<NetworKit::node> lvl, parent;
  std::vector<NetworKit::count> verticesAtLevel, prefixSum;
  NetworKit::count n, k;
  NetworKit::node l0, l1, l2;
};

struct ContractedGraph {   // "H-context"
  NetworKit::Graph H;
  PlanarGraphTools::planar_embedding_t embeddingH;
  std::vector<NetworKit::node> parentH;
  std::vector<NetworKit::count> costsH;
  NetworKit::node x;
};
```

Metody prywatne (nazwy sugerowane):

- `LevelStructure computeLevelStructure(const Graph& G, node root)` — Step 3–5
  (BFS, `verticesAtLevel`, `prefixSum`, `l1`, `k`, `l0`, `l2`), z poprawkami 1.1–1.2.
- `ContractedGraph buildContractedGraph(const Graph& G, const embedding_t& emb, const LevelStructure&)`
  — Step 6–7 (usunięcie ≥ l2, kontrakcja ≤ l0 do `x`, BFS w H, `costsH`,
  triangulacja, embedding H).
- `std::vector<node> buildInitialCycle(const ContractedGraph&, node& v1, node& w1)`
  — Step 8: pierwsza krawędź nietrzewiowa, ścieżki do LCA, budowa cyklu.
- `CycleCost computeCycleCost(const ContractedGraph&, std::vector<node>& cycle)`
  — orientacja CW + koszt inside/outside (zwraca `inside`, `outside`,
  `insideIsClockwiseArc`, zorientowany cykl). Ta sama funkcja przyda się w Step 9.

Dzięki temu `run()` staje się czytelnym szkieletem:

```cpp
void PlanarSeparator::run() {
  cleanPartitions();
  if (graph.numberOfNodes() == 0) { hasRun = true; return; }

  auto embedding  = PlanarGraphTools::findPlanarEmbedding(graph);
  auto components = NetworKit::ConnectedComponents(graph); components.run();
  auto largest    = components.extractLargestConnectedComponent(graph, false);  // raz

  if (largest.numberOfNodes() < 2 * graph.numberOfNodes() / 3) {
    findSeparatorFromComponents(components, largest);
    hasRun = true; return;
  }

  node root = firstNodeOf(largest);
  auto ls  = computeLevelStructure(largest, root);
  if (degenerate(ls)) { fallbackLevelSeparator(ls.lvl, ls.l1); hasRun = true; return; }

  auto hc  = buildContractedGraph(largest, embedding, ls);
  if (hc.H.numberOfNodes() < 3) { fallbackLevelSeparator(ls.lvl, ls.l1); hasRun = true; return; }

  node v1, w1;
  auto cycle = buildInitialCycle(hc, v1, w1);
  if (cycle.size() < 3) { fallbackLevelSeparator(ls.lvl, ls.l1); hasRun = true; return; }

  auto cc   = computeCycleCost(hc, cycle);
  // Step 9: shrinkCycleToBalanced(hc, v1, w1, cycle, cc);   // patrz osobny dokument
  auto side = markInsideOutside(hc.H, hc.embeddingH, cycle, cc.insideIsClockwiseArc);
  extractSeparatorAndPartitions(largest, ls.lvl, ls.l0, ls.l2, cycle, hc.x, side);
  hasRun = true;
}
```

Dodatkowo:

- `firstNodeOf(G)` zamiast ręcznej pętli szukającej `root`.
- W `findSeparatorFromComponents` przyjmij już wyekstrahowaną `largest`
  (usuwa drugą i trzecią ekstrakcję).
- Nazwy progów jako helpery: `costLimit() = 2 * n / 3` itp.

---

## 5. Podsumowanie priorytetów

1. **Poprawność, koniecznie:** 1.1 (`l0` break), 1.2 (`l2` fallback/`r+1`),
   1.3 (pusty graf), 1.6 (`side`), 1.4 (pokrycie całego grafu w Step 10).
2. **Robustność:** guardy przed `makeMaximalPlanar`, wejście nieplanarne,
   normalizacja kosztu dla niespójnych (1.5).
3. **Czytelność/perf:** refaktor z §4, `const auto&` zamiast kopii,
   pojedyncza ekstrakcja komponenty, `#include <numeric>`, poprawka sygnatury
   `findNumberOfVerticesAtLevel`.
