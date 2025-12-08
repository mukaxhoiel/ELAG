## ELAG — Simulation eines schiefen Kreisels mit Lagrange-Mechanik

ELAG („Euler-Lagrange Analysis for Gyroscopes“) ist ein Jupyter-Notebook, das die Bewegung eines schiefen Kreisels mithilfe der Lagrange-Mechanik symbolisch und numerisch analysiert.
Das Projekt dient als physikalische und mathematische Exploration eines Kreisels mit einem seitlich versetzten Massepunkt (Offset Gyroscope).

Ziel des Projekts:

Das Ziel ist es, die Dynamik eines Kreisels zu verstehen, der durch eine exzentrische Masse beeinflusst wird.
Dazu kombiniert das Notebook:

Symbolische Herleitungen der Lagrange-Gleichungen

Numerische Integration der Bewegungsgleichungen

Visualisierung der Winkelverläufe (ψ, φ, θ)

Darstellung von kinetischer und potenzieller Energie

Das Notebook kann als Grundlage für:

Mechanik-Studien

Simulationen von rotierenden Körpern

Lehrmaterial im Bereich klassische Mechanik

eigene physikalische Experimente
genutzt werden.

Inhalt:

Das Notebook besteht aus folgenden Abschnitten:

Symbolische Definitionen

Parameter, Winkel, Winkelgeschwindigkeiten

Position, Schwerpunkt, Drehimpuls

Kinetische & potenzielle Energie

Aufstellung der Lagrange-Funktion

Herleitung der Euler-Lagrange-Gleichungen

Numerische Lösung der Differentialgleichungen
(mit scipy.integrate.solve_ivp)

Visualisierung der Ergebnisse

Zeitverläufe der Winkel

Dynamische Interpretation

Verwendung:

Das Projekt wird direkt als Jupyter-Notebook ausgeführt.

Anforderungen:

Installiere die benötigten Python-Pakete:

pip install numpy sympy scipy matplotlib jupyter

Notebook starten:
jupyter notebook ELAG_Notebook.ipynb


Danach kannst du jeden Abschnitt einzeln ausführen, die Herleitungen ansehen und die Simulationen starten.

Abhängigkeiten

Python 3.8+

NumPy

SymPy

SciPy

Matplotlib

Jupyter Notebook

Features:

Symbolische Mechanik mit SymPy

Vollständig hergeleitete Euler-Lagrange-Gleichungen

Numerische Simulation der Dynamik des Kreisels

Visualisierung der zeitlichen Entwicklung der Winkel

Vollständig dokumentierter Rechenweg

Physikalisch motiviertes Beispiel eines realen rotierenden Systems

Lizenz — MIT License

Dieses Projekt wird unter der MIT-Lizenz veröffentlicht.
Das bedeutet:

freie Nutzung

freie Veränderung

freie Weitergabe

auch kommerzielle Nutzung erlaubt

Solange der Copyright-Hinweis erhalten bleibt.

Autor:

muka xhoiel 
GitHub: https://github.com/mukaxhoiel

## Bildrechte
Der Quellcode steht unter der MIT-Lizenz.  
Alle Bilder im Ordner `/HandNotizen` stehen jedoch unter der Lizenz 
**CC BY-NC-ND 4.0**.  
Sie dürfen heruntergeladen und privat genutzt werden, aber eine 
Weiterverwendung, Veröffentlichung oder kommerzielle Nutzung ist ohne meine 
schriftliche Zustimmung nicht gestattet.
