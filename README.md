# Molecular_dynamics_noble_gases
Molecular dynamics of one-component system with van der Waals interactions (Noble gases approximation).
Author: Krzysztof Urbanowski

ENG

Program consists of two files:
- argon.py (executable python file)
- parameters.txt (file with initital parameters)
and creates two output file:
- "result.txt"
- "avs.xyz"

Initial configuration can be changed in parameters.txt file. This file contains following variables:
* n - number of atoms along one direction (atoms creates cubic structure, so N=n^3, where N is total number of atoms)
* m - mass in atomic mass unit (u)
* a - length of unit cell in nm
* T_0 - inital temperature in K
* R - distance between the atoms (in nm) where van der Waals potential is minimum
* e - absolute value of minimal van der Waals potential (kJ/mol)
* f - elasticity coefficient of the boundary walls (kJ/(mol*nm^2))
* tau - time step in ps
* S_o - number of steps after which average values ​​are collected
* S_d - number of steps after S_o (S_d+S_o = total number of time steps)
* S_out - every S_out steps, instantaneous values are saved into result.txt
* S_xyz - every S_xyz steps, instantaneous positions of atoms are saved into avs.xyz
The configuartion given in this repository is valid for argon atoms.

After running argon.py file, user gets information about initial state in terminal (number of atomos N, T_0, kinetic energy, potential energy and total energy).
When the program is finished, user gets information about execution time and the output is stored into two files (result.txt, and avs.xyz).
* Result.txt file contains lines represting values calaculated for given time step. Each line store values in following order:
time, total energy, potential energy, temperature and pressure. First 3 values in the last line represents average values: total energy, temperature and pressure.
* avs.xyz file includes lines of x,y,z positions of atoms. Each time step is separated by number of atomos N and name of atomos in the system (in this case - Argon). This file can be opened in molecular graphics systems such as Jmol.

PL

Program składa się z dwóch plików:
- argon.py (plik wykonywalny w języku Python)
- parameters.txt (plik z początkowymi parametrami)
i tworzy dwa pliki wynikowe:
- "result.txt"
- "avs.xyz"

Początkową konfigurację można zmienić w pliku parameters.txt. Ten plik zawiera następujące zmienne:

* n - liczba atomów wzdłuż jednego kierunku (atomy tworzą strukturę kubiczną, więc N=n^3, gdzie N to całkowita liczba atomów)
* m - masa w jednostkach masy atomowej (u)
* a - długość komórki elementarnej w nm
* T_0 - początkowa temperatura w K
* R - odległość między atomami (w nm), gdzie potencjał van der Waalsa jest minimum
* e - bezwzględna wartość minimalnego potencjału van der Waalsa (kJ/mol)
* f - współczynnik elastyczności ściany granicznej (kJ/(mol*nm^2))
* tau - krok czasowy w ps
* S_o - liczba kroków, po których zbierane są wartości średnie
* S_d - liczba kroków po S_o (S_d+S_o = całkowita liczba kroków czasowych)
* S_out - co S_out kroków zapisywane są chwilowe wartości do result.txt
* S_xyz - co S_xyz kroków zapisywane są chwilowe pozycje atomów do avs.xyz
Konfiguracja zawarta w tym repozytorium jest ważna dla atomów argonu.

Po uruchomieniu pliku argon.py, użytkownik otrzymuje informacje o stanie początkowym w terminalu (liczba atomów N, T_0, energia kinetyczna, energia potencjalna i energia całkowita).
Po zakończeniu programu, użytkownik otrzymuje informacje o czasie wykonania, a wyniki są zapisywane w dwóch plikach (result.txt i avs.xyz).
* Plik result.txt zawiera linie reprezentujące obliczenia dla danego kroku czasowego. Każda linia przechowuje wartości w następującej kolejności: czas, energia całkowita, energia potencjalna, temperatura i ciśnienie. Pierwsze 3 wartości w ostatniej linii reprezentują średnie wartości: energia całkowita, temperatura i ciśnienie.
* Plik avs.xyz zawiera linie z pozycjami x, y, z atomów. Każdy krok czasowy oddzielony jest przez liczbę atomów N i nazwę atomów w systemie (w tym przypadku - Argon). Ten plik można otworzyć w systemach grafiki molekularnej, takich jak Jmol.

Numerical describtion by dr inż. Krzysztof Zberecki, Warsaw University of Technology.
Opis numeryczny autorstwa dr inż. Krzysztofa Zbereckiego, Politechnika Warszawska.

