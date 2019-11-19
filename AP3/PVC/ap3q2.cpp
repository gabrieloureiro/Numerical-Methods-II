//#EQUIPE 10
//#GABRIEL FRANÇA LOUREIRO  //#MATRÍCULA: 388835
//#ABNER DE LIMA ARAÚJO    //#MATRÍCULA: 398067
//#MN2 2019
#include "lu.hpp"
#include "PVC.hpp"
#include <iostream>

using namespace std;

//Compilação g++ 2Q.cpp LU.cpp PVC.cpp -o 2Q 
//Execução ./2Q

int main(void){

    printf ("\n>>> Entrada de valores da matrícula:\n");
    int A;
    printf ("\nA: ⇥ ");
    scanf ("%d",&A);
    int B;
    printf ("\nB: ⇥ ");
    scanf ("%d",&B);
    int C;
    printf ("\nC: ⇥ ");
    scanf ("%d",&C);
    int D;
    printf ("\nD: ⇥ ");
    scanf ("%d",&D);
    int E;
    printf ("\nE: ⇥ ");
    scanf ("%d",&E);
    int F;
    printf ("\nF: ⇥ ");
    scanf ("%d",&F);
    
    
    //Numero de Particoes
    int N = 4 + (A + B + C + D + E + F) % 4;
    //Tensao
    float T = (A + B + C);
    //Pressao
    float P = (D + E + F);

    // CCE = Condicão de Contorno Esquerda
    // CCD = Condição de Contorno Direita
    float cce = 0.2, ccd = 0.5;

    //Passo que vai ser dado do comeco da esquerda ate a direita [0,2 ate 0,5 incrementando o deltaR]
    float deltaR = (ccd - cce)/N;

    //Matriz do sistema
    MatrixXf matriz(N-1,N-1);
    matriz = PVC::getMatriz(cce,ccd,deltaR,N);

    //Vetor Y Solução
    VectorXf vetorY(N-1);

    //Vetor P
    VectorXf vetorP(N-1);
    for(int i = 0; i<N-1; i++){
        vetorP[i] = (-P/T); 
    }

    cout << endl << endl;
    cout << "N = " << N << endl;
    cout << "T = " << T << endl;
    cout << "P = " <<  P << endl;
    cout << "CCE = " << cce << endl;
    cout << "CCD = " << ccd << endl;
    cout << "∆r = " << deltaR << endl << endl;

    cout << ">>> Função dada na questão: " << endl << endl;
    cout << "T * y''(r) + T/r * y'(r) = -P " << endl;
    cout << endl << endl;
    cout << ">>> Dividindo ambos os lados por T: " << endl << endl;
    cout << "y''(r) + 1/r * y'(r) = -P/T" << endl << endl;
    printf("\nUtilizando a malha das diferenças finitas para realizar a substituição:\n\ny_i'' = (1/∆r^2)*[y_i+1 - 2_yi + y_i-1]\ny_i' = (1/2*∆r)*[y_i+1 - y_i-1]\n\n");
    printf("Reordenando a equação, conforme a malha, temos: [(1/∆r^2) - (1/2*r*∆r)]*y_i-1 + [(-2/∆r^2)]*y_i + [(1/∆r^2) + (1/2*r*∆r)]*y_i+1.\n\n");
    IOFormat CleanFmt(4, 0, ", ", "\n", "│", "│");
    cout << "⇥⇥⇥⇥⇥⇥ A * y = -P/T ⇤⇤⇤⇤⇤⇤" << endl << endl;
    cout << ">>> Matriz A: " << endl << endl <<  matriz.format(CleanFmt) << endl << endl;
    cout << ">>> Vetor y: " << endl << endl;
    cout << "y(" << cce << ") = y0 = 0" << endl;
    for(int i = 0; i<N-1; i++){cout << "y"+to_string(i+1) << endl;}
    cout << "y(" << ccd << ") = y" << N <<" = 0" << endl << endl;
    cout << ">>> Vetor P: " << endl << endl << vetorP << endl << endl;

    //Resolução do Sistema
    
    vetorY = LU::solve(matriz,vetorP,N-1);
    int posicao = 0;

    cout << "\n" <<"⇥⇥⇥⇥⇥⇥⇥⇥⇥⇥⇥⇥⇥⇥⇥ >>> AP3 QUESTÃO ➁ <<< ⇤⇤⇤⇤⇤⇤⇤⇤⇤⇤⇤⇤⇤⇤⇤:\n" << "\n" << endl;


    for(float passo=cce+deltaR; passo<ccd; passo+=deltaR){
        cout << "y(" << passo << ") = y"<< posicao+1 << " = "<< "|"<< vetorY[posicao] <<"|"<< endl;
        posicao++;
    }
    cout << endl;

    return 0;
}