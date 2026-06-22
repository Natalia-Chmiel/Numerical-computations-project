# SchrodingerSolver.jl

Pakiet w języku Julia do numerycznego rozwiązywania jednowymiarowego równania Schrödingera.

W ramach projektu zaimplementowano dwie metody numeryczne:
1. **Metodę Różnic Skończonych (FDM)** - rozwiązanie problemu własnego rzadkiej macierzy trójprzekątnej.
2. **Metodę Numerowa** - zoptymalizowany algorytm dedykowany do równań różniczkowych drugiego rzędu bez pierwszej pochodnej.

Pakiet wyznacza stany stacjonarne (dozwolone poziomy energii i funkcje falowe) dla różnych układów potencjału, testujemy je i porównujemy wyniki. Ostatecznym wynikiem projektu jest wykorzystanie tych stanów do symulacji pełnej ewolucji czasowej dla zależnego od czasu równania Schrödingera.

## Grupa
* Julia Zalewska 
* Oliwia Marut 
* Natalia Chmiel 

## Struktura projektu
* `math.ipynb` - Cała analiza matematyczna, wyprowadzenie wzorów, założenia teoretyczne oraz złożoności algorytmów.
* `codes.ipynb` - Implementacja algorytmów w języku Julia, testy wydajności, obsługa błędów oraz rozwiązanie równania Schrödingera w czasie.
* `runtests.jl` - Testy weryfikujące poprawność algorytmów.

## Uruchomienie
Wymaganymi pakietami do uruchomienia projektu są:
- `LinearAlgebra`
- `SparseArrays`
- `Plots`
- `BenchmarkTools`

