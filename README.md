# SchrodingerSolver.jl

Pakiet w języku Julia do numerycznego rozwiązywania jednowymiarowego, niezależnego od czasu równania Schrödingera. 

W ramach projektu zaimplementowano dwie metody numeryczne:
1. **Metodę Różnic Skończonych (FDM)** - rozwiązanie problemu własnego rzadkiej macierzy trójprzekątnej.
2. **Metodę Numerowa** - zoptymalizowany algorytm dedykowany do równań różniczkowych drugiego rzędu bez pierwszej pochodnej.

Następnie testujemy je na różnych potencjałach i porównujemy wyniki.

## Grupa
* Julia Zalewska 
* Oliwia Marut 
* Natalia Chmiel 

Cała analiza matematyczna, założenia i złożoności algorytmów oraz wizualizacje wyników znajdują się w pliku `math.ipynb`. W pliku `code.ipynb` znajduje się implementacja algorytmów w języku Julia.

## Uruchomienie
Wymaganymi pakietami do uruchomienia projektu są:
- `LinearAlgebra`
- `SparseArrays`
- `Plots`
- `BenchmarkTools`

