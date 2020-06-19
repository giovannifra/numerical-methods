
import math
# Esaminiamo un cuscinetto idrodinamico
P = 75000
w = 1000
rm = 0.09
l = 0.064
b = 0.046
xN = 0.038
mu = 0.04
Rt = 0.00001
# con i dati a disposizione ci troviamo la velocità della slitta e la portanza N
v = w*rm
N = P/(8*b)
print("v vale: \n", v)
print("N vale: \n", N)
# Per trovare il coefficiente caratteristico n del cuscinetto, devo ricercare gli zeri della funzione f(n).
# Siccome f(n) è non lineare, utilizzo il metodo Newton-Raphson per calcolare tali zeri.
# Prima di tutto inizializzo la variabile n arbitrariamente, per dare un valore alla prima iterazione
n = 4.5
# Definisco la tolleranza. Quando due radici di due iterazioni consecutive sono minori
# della tolleranza, il ciclo si ferma, e ci restituisce l'ultima radice trovata.
epsilon = 0.0000001
# Definisco un numero massimo di iterazioni, dopo le quali il ciclo si ferma
Niteraz = 200
# Si cercano ora gli zeri della funzione f che vale
#f = (1/2) * ((1+4*n-5*n**2+2*n*(n+2)*math.log(n))) / ((n**2-1)*math.log(n)-2*(n-1)**2) - xN/l
for i in range(Niteraz):
    f = (1/2) * ((1+4*n-5*n**2+2*n*(n+2)*math.log(n))) / ((n**2-1)*math.log(n)-2*(n-1)**2) - xN/l #funzione di cui trovare le radici
    fderiv = (((-10*n+2*(n+2)+2*n*math.log(n)+2*(n+2)*math.log(n)+4)/(2*(n**2-1)*math.log(n)-2*(n-1)**2))-((2*(n+2)*math.log(n)*n+4*n+1-5*n**2)*(2*n*math.log(n)-4*(n-1)+((n**2-1)/n))/((2*(n**2-1)*math.log(n)-2*(n-1)**2)**2)))
    nnew = n - (f/fderiv) #formula per l'iterazione
    if abs(nnew-n) < epsilon:#condizione di uscita
        n = nnew #condizione per l'iterazione
        break
    else:
        n = nnew
        continue
print("la radice vale: ", n)
print("il numero di iterazioni vale: ", i)
# Si procede col calcolo della altezza minima h2 del meato
A = 6 * ((n+1)*math.log(n)-2*(n-1)) / ((n-1)*(n**2-1))
h2 = l * (A * mu * v / N)**(1/2)
print("h2 vale: \n", h2)
# h2 è accettabile se è almeno tre volte più grande dell'altezza massima delle asperità superficiali
if h2 >= 3 * Rt:
    print("h2 accettabile")
else:
    print("h2 non valido")
# e l'altezza massima h1 del meato
h1 = n * h2
print("h1 vale: \n", h1)
# Definisco il coefficiente di attrito mediato fm
fm = ((8*mu*v)/(3*N))**(1/2) * (math.log(n)-(3*(n-1))/(2*(n+1))) / ((math.log(n)-(2*(n-1))/(n+1))**(1/2))
# fm è corretto se rientra nell'ordine del millesimo
print("fm vale: \n", fm)
if fm < 0.01:
    print("fm accettabile")
else:
    print("fm non valido")

