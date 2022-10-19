#include <stdio.h>

int globaler_wert = 0;

// Löschen Sie diese Variable nach der Bearbeitung aller Funktionen.
int TODO;

void print_pointer(int* zeiger)
{
	printf("Der Zeiger hat den Wert bzw. zeigt auf die Adresse: %p\n", 0 /* TODO.*/);
  // Hinweis: Der Compiler wird auch bei einer richtigen Antwort eine Warnung ausgeben. Fügen Sie den Type-Cast "(void *)" vor Ihrer Antwort ein, um diese zu entfernen. 
}

void print_value(int* zeiger)
{
	printf("Der Wert, auf den der Zeiger zeigt, ist: %d\n", 0 /* TODO */);
}

void set_global_value(int* zeiger)
{
  // Setzen Sie globaler_wert auf den Wert, auf den der Zeiger zeigt.
  globaler_wert = TODO;
}

void set_value(int* zeiger, int value)
{
  // Setzen Sie den Wert, auf den der Zeiger zeigt, auf den gegebenen Wert.
  TODO = value;
}

void change_pointer(int** zeiger_zeiger)
{
  // Setzen Sie den Zeiger, auf den der Zeiger zeigt, auf die Adresse von globaler_wert.
  TODO = TODO;
}

// Man kann Arrays in Funktionsköpfe verwenden:
// void change_second_value(int array_zeiger[])
// Dies ist immer identisch mit diesem Funktionskopf:
void change_second_value(int* array_zeiger)
{
  // Ändern Sie den Zeiger so, dass er auf das zweite Element des Arrays zeigt.
  array_zeiger = TODO;
  // Ändern Sie den Wert des zweiten Elements
  TODO = 200;
}

void operator_precedence()
{
  int array[4] = {1, 10, 100, 1000};
  int* zeiger_array[2] = {&array[2], &array[0]};
  
  // Entfernen Sie alle unnötigen und falschen Klammern, und fügen Sie fehlende Klammern hinzu.
  // Betrachten Sie hierfür den Begriff "Operator Precedence".
  printf("Das erste Array-Element plus 1: %d.\n", (*zeiger_array[1]) + 1);
  printf("Das zweite Array-Element: %d.\n", ((*zeiger_array)[1]) + 1);
  printf("Das vierte Array-Element plus 1: %d.\n", *((zeiger_array[1]) + 1));
}

int main (void)
{
  int lokale_werte[2] = {0, 100};
  int* zeiger = &lokale_werte[0];
  *zeiger = 5;

  // Den Zeiger lesend verwenden.
  print_pointer(zeiger);
  print_value(zeiger);
  printf("\n");
  
  // Einen Zeiger in einer Zuweisung verwenden.
  printf("Der globale Wert (aktuell: %d) wird auf den ersten lokalen Wert (%d) gesetzt.\n", globaler_wert, lokale_werte[0]);
  set_global_value(zeiger);
  printf("Der globale Wert ist jetzt: %d.\n", globaler_wert);
  printf("\n");
  
  // Einen Wert mit einem Zeiger ändern
  printf("Der erste lokale Wert wird auf 10 gesetzt.\n");
  set_value(zeiger, 10);
  print_pointer(zeiger);
  print_value(zeiger);
  printf("\n");
  
  // Den Zeiger ändern
  printf("Der Zeiger soll nun auf den globalen Wert zeigen.\n");
  printf("Der globale Wert ist momentan %d.\n", globaler_wert);
  printf("Die Addresse vom globalen Wert ist momentan %p.\n", (void *) &globaler_wert);
  change_pointer(&zeiger);
  print_pointer(zeiger);
  print_value(zeiger);
  printf("Jetzt wird mit dem Zeiger der globale Wert auf 20 geändert.\n");
  set_value(zeiger, 20);
  printf("Globaler Wert:\n");
  print_pointer(zeiger);
  print_value(zeiger);
  printf("Lokaler Wert:\n");
  print_pointer(&lokale_werte[0]);
  print_value(&lokale_werte[0]);
  printf("\n");
  
  // Zeiger als Array-Zeiger
  printf("Arrays sind in gewisser Weise Zeiger auf Werte.\n");
  printf("Der zweite Wert des Arrays ist gerade %d.\n", lokale_werte[1]);
  printf("Der zweite Wert des Arrays wird auf 200 geändert.\n");
  change_second_value(lokale_werte);
  printf("Jetzt ist der zweite Wert %d.\n", lokale_werte[1]);
  printf("\n");
  
  // Operator Precedence
  printf("Schließlich ein Exkurs in die Operator Precedence:\n");
  operator_precedence();
  
  return 0;
}
