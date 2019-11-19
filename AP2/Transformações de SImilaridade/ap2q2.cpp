//#EQUIPE 10
//#GABRIEL FRANÇA LOUREIRO  //#MATRÍCULA: 388835
//#ABNER DE LIMA ARAÚJO    //#MATRÍCULA: 398067
//#MN2 2019
#include "lib.hpp"

int main(){

	printf ("\n>>> Insira a matrícula a ser utilizada:\n");
	int a;
	printf ("\nA: ⇥ ");
  	scanf ("%d",&a);
	int b;
	printf ("\nB: ⇥ ");
  	scanf ("%d",&b);
	int c;
	printf ("\nC: ⇥ ");
  	scanf ("%d",&c);
	int d;
	printf ("\nD: ⇥ ");
  	scanf ("%d",&d);
	int e;
	printf ("\nE: ⇥ ");
  	scanf ("%d",&e);
	int f;
	printf ("\nF: ⇥ ");
  	scanf ("%d",&f);




	MatrixXd A(5,5);
	A <<
			95, 10, 0, 5, 2,
			10, 45, 10, 7, 4,
			0, 10, 30, 9, 1,
			5, 7, 9, 20, 0.5,
			2, 4, 1, 0.5, 80;

	//Erro
	double E = 0.000001;

	int rows = A.rows();
	int cols = A.cols();
	
	//--------------------------- Achar Matriz U pelo Metodo de Jacobi ------------------------------------
	MatrixXd A_U(rows, rows);
	MatrixXd H_U(rows, rows);
	MatrixXd It_U(rows, rows);
	MatrixXd valueU(rows, rows);
	MatrixXd vectorU(rows, rows);

	A_U = A*A.transpose();
	tie(It_U,H_U) = HouseHolder(A_U);
	tie(valueU, vectorU) = jacobi(It_U, E, H_U);

	around(It_U);
	around(H_U);
	around(valueU);
	around(vectorU);

	ordenar(valueU, vectorU);


	IOFormat CleanFmt(4, 0, ", ", "\n", "│", "│");
	cout << "\n" <<"######################## >>> Método de Jacobi <<< #######################" << "\n" << endl;
	cout << ">>> Matriz A*A^T: " << endl << A_U.format(CleanFmt) << "\n" << endl;
	cout << ">>> Matriz Tridiagonal: " << endl << It_U.format(CleanFmt) << "\n" << endl;
	cout << ">>> Matriz de HouseHolder Acumulada: " << endl << H_U.format(CleanFmt) << "\n" << endl;
	cout << ">>> Matriz de Autovalores: "<< endl << valueU.format(CleanFmt) << "\n" << endl;
	cout << ">>> Matriz de Autovetores: " << endl << vectorU.format(CleanFmt) << "\n" << endl;

	

	//------------------------------- Achar Matriz V pelo Metodo QR ------------------------------------
	MatrixXd A_V(cols, cols);
	MatrixXd H_V(cols, cols);
	MatrixXd It_V(cols, cols);
	MatrixXd valueV(cols, cols);
	MatrixXd vectorV(cols, cols);
	
	A_V = A.transpose()*A;
	tie(It_V,H_V) = HouseHolder(A_V);
	tie(valueV, vectorV) = QR(It_V, E, H_V);
	
	around(It_V);
	around(H_V);
	around(valueV);
	around(vectorV);

	ordenar(valueV, vectorV);


	cout << "\n" <<"######################## >>> Método QR <<< #######################" << "\n" << endl;
	cout << ">>> Matriz A^T*A: " << endl << A_V.format(CleanFmt) << "\n" << endl;
	cout << ">>> Matriz Tridiagonal: " << endl << It_V.format(CleanFmt) << "\n" << endl;
	cout << ">>> Matriz de HouseHolder Acumulada: " << endl << H_V.format(CleanFmt) << "\n" << endl;
	cout << ">>> Matriz de Autovalores: "<< endl << valueV.format(CleanFmt) << "\n" << endl;
	cout << ">>> Matriz de Autovetores: " << endl << vectorV.format(CleanFmt) << "\n" << endl;



	//------------------------------------------- RESULTADO ---------------------------------------
	

	// Sigma
	MatrixXd Sigma(rows, cols);
	Sigma = MatrixXd::Zero(rows,cols);
	for (int i = 0; i < rows; i++) 
		Sigma(i,i) = valueU(i,i);
    sqrt_diagonal(Sigma);

    consertaSinal(A, vectorV, vectorU, Sigma);

    // U * Σ * Vt
    MatrixXd aux(rows, cols);
    aux = vectorU*Sigma*vectorV.transpose();
    around(aux);


	cout << "\n" <<"⇥⇥⇥⇥⇥⇥⇥⇥⇥⇥⇥⇥⇥⇥⇥ >>> AP2 QUESTÃO ➁ <<< ⇤⇤⇤⇤⇤⇤⇤⇤⇤⇤⇤⇤⇤⇤⇤:\n" << "\n" << endl;
	printf("Definição: A = U*Σ*(V^t)\n\n>>> Prove que: [A*(A^T)]*U = U*Σ*(Σ^t)\n\nPela definição, temos que A*(A^t) = U*Σ*(V^t)*V*(Σ^t)*(U^t), mas, como V é ortogonal((V^t)*V = I), temos A*(A^t) =  U*Σ*(Σ^t)*(U^t).\n\nMultiplicando ambos os lados por U, ficamos com [A*(A^t)]*U = U*Σ*(Σ^t)*(U^t)*U, mas, como U é ortogonal((U^t)*U = I), temos [A*(A^t)]*U = U*Σ*(Σ^t).\n\n");
	printf("Definição: A = U*Σ*(V^t)\n\n>>> Prove que: [(A^T)*A]*V = V*(Σ^t)*Σ\n\nPela definição, temos que (A^t)*A = V*(Σ^t)*(U^t)*U*Σ*(V^t), mas, como U é ortogonal((U^t)*U = I), temos (A^t)*A =  V*(Σ^t)*Σ*(V^t).\n\nMultiplicando ambos os lados por V, ficamos com [(A^t)*A]*V = V*(Σ^t)*Σ*(V^t)*V, mas, como V é ortogonal((V^t)*V = I), temos [(A^t)*A]*V = V*(Σ^t)*Σ.\n\n");
	cout << ">>> Matriz U: " << endl << vectorU.format(CleanFmt) << "\n" << endl;
	cout << ">>> Matriz Σ: " << endl << Sigma.format(CleanFmt) << "\n" << endl;
	cout << ">>> Matriz V: " << endl << vectorV.format(CleanFmt) << "\n" << endl;
	cout << ">>> Matriz A: " << endl << A.format(CleanFmt) << "\n" << endl;
	cout << ">>> Matriz U*Σ*(V^t): " << endl << aux.format(CleanFmt) << "\n" << endl;


	return 0;
}