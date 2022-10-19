#include <stdio.h>

// Definieren Sie ein 3x3-Array Namens map, das Werte vom Typ double enthält
double map[3][3];

// Die Funktion set_temperature soll an Position [x, y] den Wert dir in das Array map eintragen
// Überprüfen Sie x und y, um mögliche Arrayüberläufe zu verhindern
void set_temperature(int x, int y, double temperature)
{
  if (x < 0 || x >= 3 || y < 0 || y >= 3)
    return;
  map[x][y] = temperature;
}

/**
 * Gibt die Temperatur an der [x,y]-Koordinate aus, oder 0, wenn sie außerhalb des Arrays liegt.
 */
double get_temperature(int x, int y)
{
  if (x < 0 || x >= 3 || y < 0 || y >= 3)
    return 0.0;
  return map[x][y];
}

// Die Funktion show_map soll das Array in Form einer 3x3-Matrix ausgeben
void show_map(void)
{
  for (int y = 0; y < 3; y++)
  {
    for (int x = 0; x < 3; x++)
    {
      printf("%f\t", map[x][y]);
    }
    printf("\n");
  }
}

/**
 * Summiert alle 9 Werte um [x_center,y_center]
 * Alle Werte außerhalb des Arrays sind 0.
 */
double sum_surrounding(int x_center, int y_center)
{
  double sum = 0;
  for (int y = y_center - 1; y <= y_center + 1; y++)
  {
    for (int x = x_center - 1; x <= x_center + 1; x++)
    {
      sum += get_temperature(x, y);
    }
  }
  return sum;
}

// Die Funktion set_to_average soll an Position [x, y] den Durchschnitt der 8 umgebenen
// Temperaturen in das Array map eintragen.
// Für Werte außerhalb des Arrays soll der Wert 0 angenommen werden.
// Verwenden Sie hierfür auch die Funktion set_temperature.
/**
 * Trägt den Durschnitt der 8 Werte um [x,y] bei [x,y] ein.
 */
void set_to_average(int x, int y)
{
  // [x,y] auf 0 setzen, damit der Wert nicht zur Summe beiträgt.
  set_temperature(x, y, 0);
  // Summe berechnen
  double sum = sum_surrounding(x, y);
  // Durchschnittswert eintragen
  set_temperature(x, y, sum / 8.0);
}

// In dieser Funktion darf nichts verändert werden!
int main (void)
{
	set_temperature(0, 1, 40);
	set_temperature(1, 0, 160);
	set_temperature(1, 4, 75);
	set_temperature(1, 2, 80);
	set_temperature(2, 1, 120);

	show_map();

	set_temperature(0, 0, 20.5);
	set_temperature(0, 2, 14.8);
	set_temperature(0, 2, 22.7);
	set_temperature(2, 0, 100.2);
	set_temperature(2, 2, 20.6);
	set_temperature(2, 2, 200.2);
	set_temperature(1, 3, 201.06);
	set_temperature(1, 1, 50.5);

	show_map();
  
  set_to_average(0,0);
  set_to_average(2,0);
  set_to_average(1,2);
  
  show_map();

	return 0;
}
