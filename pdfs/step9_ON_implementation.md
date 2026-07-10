# Step 9 separatora planarnego (Lipton–Tarjan) — implementacja w O(n)

Dokument opisuje, jak zaimplementować **Step 9** algorytmu z pracy
_R. J. Lipton, R. E. Tarjan, „A Separator Theorem for Planar Graphs", SIAM J.
Appl. Math. 36(2), 1979_ w prawdziwym czasie **O(n)** (a nie O(n²)), dopasowaną
do stanu kodu w `cpp/separator/PlanarSeparator.cpp` i `include/separator/PlanarSeparator.hpp`.

Klucz do O(n) jest dokładnie taki, jak w akapicie _„Proof of Step 9 time bound"_:
**nigdy nie odbudowujemy cyklu ani nie liczymy kosztu od zera** — koszt „inside"
utrzymujemy przyrostowo, a w przypadku B używamy naprzemiennego skanu, który
kosztuje tyle, ile _mniejsza_ strona.

> Sygnatury `applyStep9CaseB` / `markInsideOutside` w nagłówku są skrojone pod
> wersję wolną. Dla O(n) stan to `(v1, w1, curDart, insideCost, onCycle[])`, a
> nie `cycle` + pełne `side` w każdej iteracji — nagłówek trzeba odpowiednio
> zmienić.

---

## 1. Klucz do O(n) — wzory na koszt (bez przeliczania)

Niech `total` = liczba wierzchołków w `H` (koszt jednostkowy 1; ewentualnie
`total = costsH[x] - 1`, bo wierzchołek `x` ma mieć koszt 0). `insideCost`
trzymamy jako jedną liczbę i aktualizujemy przyrostowo.

### Przypadek A — jedna z krawędzi (v1, y), (y, w1) jest krawędzią drzewa

Nowy cykl = stary minus trójkąt (v1, y, w1); jedyna zmiana to `y` wchodzi na cykl:

```
insideCost_new = insideCost_old - cost(y)
```

Wszystko O(1). (Lemat 2, przypadek 3b.)

### Przypadek B — obie krawędzie (v1, y), (y, w1) są nietrzewiowe

Idziemy po `parentH` od `y`, aż trafimy na wierzchołek cyklu `z`. Ścieżka
`y -> z` jest wspólną granicą obu podcykli, więc:

```
inside_1 + inside_2 = insideCost_old - P
P = suma cost(u) dla u na ścieżce (y -> z), u != z      // wliczając y
```

Naprzemiennym skanem liczymy koszt **mniejszej** strony (np. `inside_1`), a drugą
dostajemy odejmowaniem:

```
inside_2       = insideCost_old - P - inside_1
insideCost_new = max(inside_1, inside_2)
```

Nowy `(v1, w1)` = krawędź nietrzewiowa strony o większym koszcie. (Lemat 2,
przypadek 3c.)

---

## 2. Preprocessing O(n) — struktura półkrawędzi (apex w O(1))

Bez tego znajdowanie trójkąta byłoby O(deg) i suma nie wyszłaby O(n).
`embeddingH` (rotacje CW) triangulowanego `H` zamieniamy na nawigację po ścianach:

```cpp
NetworKit::count nb = H.upperNodeIdBound();
std::vector<std::unordered_map<NetworKit::node,int>> idxOf(nb);
H.forNodes([&](NetworKit::node u){
    auto &rot = embeddingH[u];
    for (int i = 0; i < (int)rot.size(); ++i) idxOf[u][rot[i]] = i;
});
auto deg  = [&](NetworKit::node u){ return (int)embeddingH[u].size(); };
using Dart = std::pair<NetworKit::node,int>;              // u -> embeddingH[u][idx]
auto head = [&](Dart d){ return embeddingH[d.first][d.second]; };
auto twin = [&](Dart d)->Dart{ NetworKit::node v=head(d); return {v, idxOf[v][d.first]}; };
// następna półkrawędź wzdłuż ściany po USTALONEJ stronie.
// UWAGA: kierunek (+1 vs -1) zależy od orientacji embeddingu — jeśli obchodzi
// ścianę zewnętrzną, zmień na (t.second - 1 + deg)%deg.
auto nextInFace = [&](Dart d)->Dart{ Dart t=twin(d); return {t.first,(t.second+1)%deg(t.first)}; };
auto apexOf = [&](Dart d){ return head(nextInFace(d)); };  // trzeci wierzchołek trójkąta po stronie d
```

---

## 3. Stan utrzymywany między iteracjami

```cpp
std::vector<char> onCycle(nb, 0);
for (auto v : cycle) onCycle[v] = 1;                 // z Step 8

int total = static_cast<int>(H.numberOfNodes());     // lub costsH[x]-1
int insideCost = inside;                             // z Step 8

// curDart = półkrawędź krawędzi nietrzewiowej (v1,w1), której LEWA ściana leży WEWNĄTRZ.
Dart curDart = {v1, idxOf[v1][w1]};
if (onCycle[ apexOf(curDart) ] /*apex musi być ściśle wewnątrz*/) curDart = twin(curDart);
// (praktycznie: wybierz z {v1->w1, w1->v1} ten, którego apex nie jest na cyklu
//  i leży po stronie 'insideIsClockwiseArc' wyliczonej w Step 8)
auto vertexCost = [&](NetworKit::node){ return 1; };  // koszt jednostkowy
```

---

## 4. Pętla Step 9

```cpp
while (insideCost * 3 > 2 * total) {
    NetworKit::node y = apexOf(curDart);
    NetworKit::node a = curDart.first, b = head(curDart);   // bieżąca krawędź nietrzewiowa (a,b) = (v1,w1)

    bool tAy = (parentH[a] == y || parentH[y] == a);   // (a,y) drzewowa?
    bool tyB = (parentH[b] == y || parentH[y] == b);   // (y,b) drzewowa?

    if (tAy || tyB) {
        // ---- Przypadek A: O(1) ----
        onCycle[y] = 1;
        insideCost -= vertexCost(y);
        // trójkąt = LEWA ściana curDart: d0=curDart(a->b), d1=nextInFace(d0)(b->y), d2=(y->a)
        Dart d1 = nextInFace(curDart);        // b -> y
        Dart d2 = nextInFace(d1);             // y -> a
        // nowa krawędź nietrzewiowa = ta NIE-drzewowa; nowe wnętrze jest po DRUGIEJ
        // stronie tej krawędzi niż trójkąt, więc bierzemy twin.
        curDart = tAy ? twin(d1)   // (a,y) drzewowa -> nowa nietrzewiowa (y,b): dart y->b
                      : twin(d2);  // (y,b) drzewowa -> nowa nietrzewiowa (a,y): dart a->y
    } else {
        // ---- Przypadek B ----
        applyStep9CaseB_On(/* aktualizuje a,b,curDart,insideCost,onCycle */);
    }
}
```

---

## 5. Przypadek B w O(mniejsza strona + |path|)

Krok po kroku, dokładnie jak w pracy.

### 5.1. Ścieżka y -> z i koszt P

```cpp
int P = 0;
std::vector<NetworKit::node> zpath;
NetworKit::node z = y;
while (!onCycle[z]) { zpath.push_back(z); P += vertexCost(z); z = parentH[z]; }
for (auto u : zpath) onCycle[u] = 1;          // teraz granice OBU podcykli są onCycle
```

Po zamarkowaniu ścieżki `y -> z` wszystkie wierzchołki brzegowe **obu** podcykli
są `onCycle` (dowód: L1 = LCA(v1, y) to albo `z`, albo stare LCA `L`, a ramiona
podcykli są podścieżkami starego cyklu lub ścieżki `y -> z`).

### 5.2. Darty startowe podcykli (wnętrze po lewej)

```cpp
Dart d1 = nextInFace(curDart);   // b -> y
Dart d2 = nextInFace(d1);        // y -> a
Dart start1 = twin(d2);          // krawędź (a,y): dart a->y, lewa ściana = wnętrze podcyklu 1
Dart start2 = twin(d1);          // krawędź (y,b): dart y->b, lewa ściana = wnętrze podcyklu 2
```

### 5.3. Naprzemienny skan ścian

Granica podcyklu = krawędź nietrzewiowa (nie przekraczamy jej twinem na zewnątrz)

- krawędzie drzewowe o **obu** końcach `onCycle`.

```cpp
auto isBoundary = [&](Dart d, Dart nontreeDart){
    if (d == nontreeDart || twin(d) == nontreeDart) return true;
    NetworKit::node u = d.first, v = head(d);
    bool tree = (parentH[u]==v || parentH[v]==u);
    return tree && onCycle[u] && onCycle[v];
};
```

Skanery uruchamiamy w pętli po jednej ścianie z każdej strony, aż któryś skończy
kolejkę:

```cpp
Scan s1(start1, /*nontree=*/twin(start1));   // wnętrze podcyklu (a,y)
Scan s2(start2, /*nontree=*/twin(start2));   // wnętrze podcyklu (y,b)
bool firstDone1;
while (true) {
    s1.stepOneFace();   // odwiedza 1 ścianę, dolicza NOWE wierzchołki wnętrza
    s2.stepOneFace();
    if (s1.finished()) { firstDone1 = true;  break; }
    if (s2.finished()) { firstDone1 = false; break; }
}
int insideSmall = firstDone1 ? s1.cost : s2.cost;
int inside1 = firstDone1 ? insideSmall : (insideCost - P - insideSmall);
int inside2 = firstDone1 ? (insideCost - P - insideSmall) : insideSmall;
```

`Scan::stepOneFace` to BFS w dualu: z bieżącej ściany dla każdej z 3 półkrawędzi
`d` sprawdź `!isBoundary(d, nontree)` → przejdź `twin(d)` do sąsiedniej ściany
(jeśli nieodwiedzona); dla każdej nowej ściany dolicz `vertexCost` jej
wierzchołków, które nie są `onCycle` i nie były jeszcze zliczone (znacznik + lista
do resetu po zakończeniu). Koszt jednego skanu = O(#ścian wnętrza), a
naprzemienność ogranicza go do mniejszej strony.

### 5.4. Wybór strony i sprzątanie

```cpp
bool takeSide1 = (inside1 >= inside2);
insideCost = takeSide1 ? inside1 : inside2;
curDart    = takeSide1 ? start1 : start2;   // dart nietrzewiowej wybranej strony, wnętrze po lewej

// odmarkuj ramię strony ODRZUCONEJ (od jej endpointu w górę aż do z, wyłącznie):
NetworKit::node drop = takeSide1 ? b /*stare w1*/ : a /*stare v1*/;
for (NetworKit::node u = drop; u != z && onCycle[u]; u = parentH[u]) onCycle[u] = 0;
```

---

## 6. Dlaczego to jest O(n) (mapowanie na dowód z pracy)

- Przypadek A: O(1).
- Przypadek B: O(|path(y -> z)|) + O(min(wnętrze_1, wnętrze_2)) (skan)
  - O(|odrzucone ramię|) (odmarkowanie).
- Każdy wierzchołek ścieżki `y -> z` (poza `z`) jest wewnątrz bieżącego, ale poza
  wszystkimi kolejnymi cyklami → policzony raz.
- Każda zeskanowana ściana/wierzchołek mniejszej strony trafia poza wszystkie
  kolejne cykle → policzone raz. Odrzucone ramię również nigdy nie wraca na cykl.
- Stąd suma po wszystkich iteracjach = O(n). Dokładnie akapit _„For every two
  edges scanned, at least one edge is inside the current cycle but outside all
  subsequent cycles."_

---

## 7. Po pętli — side[] raz, potem Step 10

`markInsideOutside` (flood fill od wnętrza, ściany = `onCycle`) wołasz **raz**,
O(n) — nie psuje bilansu:

```cpp
std::vector<Side> side = markInsideOutside(H, embeddingH, /*cycle z onCycle*/, insideCW);
extractSeparatorAndPartitions(G, lvl, l0, l2, /*cycle*/, x, side);
```

---

## 8. Do zweryfikowania przy wklejaniu

Inaczej dostaniesz ścianę zewnętrzną albo złą stronę:

1. Kierunek `nextInFace` (`+1` vs `-1`) względem orientacji Twojego `embeddingH`.
2. Wybór `curDart` po Step 8 tak, by lewa ściana faktycznie była wnętrzem
   (`insideIsClockwiseArc`).
3. W przypadku B poprawne zapamiętanie starych `v1`, `w1` do odmarkowania ramienia.
