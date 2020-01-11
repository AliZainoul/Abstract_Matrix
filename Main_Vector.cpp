#include <iostream>
#include "Vector.hpp"
using namespace std;

// La fonction fill() est sujette au "copy elision"
Vector fill()
{
	cout << "Entrée fill()"<<"\n";
	Vector vec={8,8,8,8};
	cout << "Sortie fill()"<<"\n";
	return vec;
}
// La fonction change() n'est pas sujette au "copy elision"
//  (il est interdit au compilateur
//  d'éluder la copie d'un argument passé par valeur
//  vers un objet temporaire en return)
Vector change(Vector tmp)
{
	cout << "Entrée change()"<<"\n";
	tmp[3]=500;
	cout << "Sortie change()"<<"\n";
	return tmp;
}

int main()
{
  Vector v3(4,1);
  Vector v4={-4.3,4.5,-4.6,-4.7};
  Vector v5(5,5);

  // Initialisation ==> Copy Constructor
  cout<< "Test du constructeur de copie 1"<<"\n";
  Vector w1{v3};
  cout << w1;
  cout << "La norme L2 du vecteur: "<< w1.twonorm()<<"\n";
  cout<< "                        "<<"\n";

  // Affectation ==> Copy Assignment
  cout<< "test assign assignment 1b"<<"\n";
  Vector w1b(4);
  w1b = v4;
  cout << w1b;
	cout << "La norme Linfinity du Vecteur: "<< w1b.maxnorm()<<"\n";
  cout<< "                        "<<"\n";

  // Initialisation ==> Copy constructor
  cout<< "Test du constructeur de copie 2"<<"\n";
  Vector w2(v3);
  cout << w2;
  cout<< "                        "<<"\n";


  // Initialisation
  // le constructeur de copy n'est pas appelé
  // c'est une initialisation !
	cout<< "Test du constructeur de copie 3"<<"\n";
  Vector v6=Vector(6,6);
  cout << v6;
  cout<< "copy elision !"<<"\n";
  cout<< "                        "<<"\n";

	// le constructeur de copie n'est pas appelé car le vecteur issu de fill ne sera jamais utilisé "seul"
	cout<< "Test du constructeur de copie 4"<<"\n";
	Vector v7=fill();
	cout << v7;
	cout<< "copy elision !"<<"\n";
	cout<< "                        "<<"\n";


  // Affectation ==> move assignment car Vector(6,6) ne sera jamais utilisé "seul"
  cout<< "test copy/move constructeur 3b"<<"\n";
  Vector v6b(6);
  v6b = Vector(6,6);
  cout << v6b;
  cout<< "                        "<<"\n";

  // cout<< "test copy constructeur 5"<<"\n";
  // Vector w3=v5;
  // cout << w3;
  // cout<< "pas d'appel de l'opÃ©rateur d'assignement car initialisation!"<<"\n";
  // cout<< "                        "<<"\n";

  cout<< "test move constructor"<<"\n";
  Vector v8(8,8);
  Vector v8b=change(v8);
  cout << v8b;
  cout<< "maintenant, move constructor est appelÃ© automatiquement"<<"\n";
  cout<< "                        "<<"\n";

  cout<< "test move assignment"<<"\n";
  Vector v9(9,9);
  cout << v9;
  v9 = fill();
  cout << v9;

  cout<< "move assignment est appelé automatiquement"<<"\n";
  cout<< "car il s'agit d'une affectation et non initialisation"<<"\n";
  cout<< "                        "<<"\n";


}
