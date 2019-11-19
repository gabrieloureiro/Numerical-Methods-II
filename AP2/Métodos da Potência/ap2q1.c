//#EQUIPE 10
//#GABRIEL FRANÇA LOUREIRO  //#MATRÍCULA: 388835
//#ABNER DE LIMA ARAÚJO    //#MATRÍCULA: 398067
//#MN2 2019
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//Definição das structs
struct vetorAutovalores{
    double   autoValor;
    double *autoVetor;
};

struct ConjvetorAutovalores{
    double  **autoValores;
    double **autoVetores;
};

//Definição do tamanho da matriz e do vetor

int colunas    = 5;
int linhas      = 5;
int tamanhoVetor = 5;


/// Função para alocação do espaço de um vetor
/// \return vetor com o espaço alocado
double* alocarVetor(){
    double *vetor = (double*) malloc(tamanhoVetor * sizeof(double));
    return vetor;
}

/// Função para inicializar o vetor com 0
/// \param vetor
void inicializarVetor(double *vetor) {
    for (int l = 0; l < linhas ; ++l) {
        vetor[l] = 0;
    }
}

/// Função para alocação do espaço para uma Matriz
/// \param linhas = quantidade de linhas
/// \param colunas = quantidade de colunas
/// \return a matriz com o espaço alocado
double** alocarMatriz(int linhas, int colunas){
    double **matriz = (double**)malloc(linhas * sizeof(double*));

    for (int i = 0; i < linhas; i++){
        matriz[i] = (double*) malloc(colunas * sizeof(double));
    }
    return matriz;
}

/// Imprime a Matriz
/// \param matriz
void printarMatriz(double **matriz){
    for (int l = 0; l < linhas; ++l) {
        for (int m = 0; m < colunas; ++m) {
            printf("%f", matriz[l][m]);
            printf("\t");
        }
        printf("\n");
    }
    printf("\n");
}

/// Imprime Vetor
/// \param vetor
void printVetor(double *vetor){
    for (int i = 0; i < linhas; ++i) {
        printf("%f", vetor[i]);
        printf("\n");
    }
    printf("\n");
}

/// Normalização do vetor
/// \param vetor
/// \return a norma do vetor
double normalizacaoEuclidiana(double *vetor) {
    double normalizacao = 0;
    for (int j = 0; j < tamanhoVetor; ++j) {
        normalizacao += pow(vetor[j], 2);
    }
    normalizacao = sqrt(normalizacao);
    return normalizacao;
}

/// Faz a normalização euclidiana de uma matriz
/// \param matriz
/// \return a norma da matriz
double normalizacaoMatrizEuclidiana(double **matriz){
    double normalizacao = 0;
    for (int i = 0; i < linhas; ++i) {
        for (int j = 0; j < colunas; ++j) {
            normalizacao += pow(matriz[i][j], 2);
        }
    }
    normalizacao = sqrt(normalizacao);
    return normalizacao;
}

/// Multiplica matriz e vetor
/// \param matriz
/// \param vetor
/// \return o vetor resultante da multiplicação
double* matriz_x_vetor(double **matriz, const double *vetor){
    double *vetorResultante = (double*) malloc(tamanhoVetor * sizeof(double));
    for (int i = 0; i < linhas ; ++i) {
        vetorResultante[i] = 0;
        for (int j = 0; j < colunas; ++j) {
            vetorResultante[i] += vetor[j] * matriz[i][j];
        }
    }
    return vetorResultante;
}

/// Faz a matriz diagonal de uma dada matriz
/// \param matriz
/// \return matriz diagonal
double** criarmatrizDiagonal(double **matriz){
    double **matriz_diagonal = alocarMatriz(linhas, colunas);
    for (int l = 0; l < linhas ; ++l) {
        for (int m = 0; m < colunas; ++m) {
            if (l == m){
                matriz_diagonal[l][m] = matriz[l][m];
            }
            else{
                matriz_diagonal[l][m] = 0.0;
            }
        }
    }
    return matriz_diagonal;
}

/// Normaliza o vetor
/// \param vetor
/// \return o vetor normalizado
double* normalizarVetor(double *vetor) {
    double *normVetor = (double*) malloc(tamanhoVetor * sizeof(double));
    double     normalizacao = normalizacaoEuclidiana(vetor);
    for (int i = 0; i < tamanhoVetor; ++i) {
        normVetor[i] = vetor[i]/normalizacao;
    }
    return normVetor;
}

/// Multiplicação de vetores
/// \param vetor1
/// \param vetor2
/// \return o resultado da multiplicacao
double multiplicarVetores(const double *vetor1, const double *vetor2) {
    double result = 0;
    for (int i = 0; i < tamanhoVetor ; ++i) {
        result += vetor1[i] * vetor2[i];
    }
    return result;
}

/// Multiplicacao de um vetor transposto por um vetor normal
/// \param vetor1
/// \param vetor2
/// \return Matriz resultante
double** vetorTransposto_x_vetorNormal(const double *vetor1, const double *vetor2){
    double **matrizResultante  = alocarMatriz(linhas, colunas);
    for (int i = 0; i < colunas; ++i) {
        for (int j = 0; j < linhas ; ++j) {
            matrizResultante[j][i] = vetor1[j] * vetor2[i];
        }
    }
    return matrizResultante;
}

/// Subtração entre dois vetores
/// \param vetor1
/// \param vetor2
/// \return Vetor resultante
double* subtrairVetor(const double *vetor1, const double *vetor2){
    double *vetorResultante = alocarVetor();
    for (int i = 0; i < linhas; ++i) {
        vetorResultante[i] = vetor1[i] - vetor2[i];
    }
    return vetorResultante;
}

/// Transforma uma matriz numa matriz identidade
/// \param matriz
void matrizIdenditade(double **matriz) {
    for (int l = 0; l < linhas ; ++l) {
        for (int m = 0; m < colunas; ++m) {
            if (l == m){
                matriz[l][m] = 1.0;
            }
            else{
                matriz[l][m] = 0.0;
            }
        }
    }
}

/// Função para calcular a matriz resultante da multiplicação de uma matriz por um valor escalar
/// \param matriz
/// \param value
/// \return Matriz resultado da multiplicação
double** matriz_x_escalar(double **matriz, double value){
    double **matrizResultante = alocarMatriz(linhas, colunas);
    for (int l = 0; l < linhas ; ++l) {
        for (int m = 0; m < colunas; ++m) {
            matrizResultante[l][m] = matriz[l][m] * value;

        }
    }
    return matrizResultante;
}

/// Função para multiplicar matrizes
/// \param matriz1
/// \param matriz2
/// \return Matriz resultado da multiplicação
double** multiplicarMatriz(double **matriz1, double **matriz2){
    double **matrizResultante = alocarMatriz(linhas, colunas);
    for (int i = 0; i < linhas ; ++i) {
        for (int j = 0; j < colunas ; ++j) {
            matrizResultante[i][j] = 0;
            for (int k = 0; k < colunas; ++k) {
                matrizResultante[i][j] += matriz1[i][k] * matriz2[k][j];
            }
        }
    }
    return matrizResultante;
}

/// Transpõe uma matriz
/// \param matriz
/// \return a matriz transposta
double** matrizTransposta(double **matriz){
    double **matrizTransposta = alocarMatriz(linhas, colunas);
    for(int i=0; i < linhas; ++i)
        for(int j=0; j < colunas; ++j) {
            matrizTransposta[j][i] = matriz[i][j];
        }
    return matrizTransposta;
}

/// Função que calcula a matriz resultante da subtração de duas matrizes
/// \param matriz1
/// \param matriz2
/// \return Matriz resultado da subtração
double** subtrairMatriz(double **matriz1, double **matriz2){
    double **matrizResultante = alocarMatriz(linhas, colunas);
    for (int l = 0; l < linhas ; ++l) {
        for (int m = 0; m < colunas; ++m) {
            matrizResultante[l][m] = matriz1[l][m] - matriz2[l][m];

        }
    }
    return matrizResultante;
}

/// Função para checar os vetores e identificar se é necessário mudar o sinal
/// \param indice
/// \param vetorP
/// \param vetor_aux
/// \param vetorPnormalizacao
void verificarSinal(int indice, const double *vetorP, double *vetorP_aux, double vetorP_normalizado) {
    if(vetorP[indice + 1] > 0){
        vetorP_aux[indice + 1] = vetorP_normalizado * -1.0;
    }

    else {
        vetorP_aux[indice + 1] = vetorP_normalizado;
    }
}

/// Calcula o maior autovalor pelo metodo da potencia regular
/// \param matriz        = Matriz Original
/// \param vetorInicial = Vetor chute
/// \param tolerancia     = Tolerancia para o erro
/// \return struct com o maior autovalor e o autovetor correspondente
struct vetorAutovalores potenciaRegular(double **matriz, double *vetorInicial, double tolerancia){
    struct vetorAutovalores resultado;
    resultado.autoVetor = alocarVetor();
    double *q = normalizarVetor(vetorInicial);
    double erro;
    double autoValor;
    double autoValor_aux = 0;
    double *x = matriz_x_vetor(matriz, q);

    do{
        q = normalizarVetor(x);
        x = matriz_x_vetor(matriz, q);
        autoValor = multiplicarVetores(q, x);
        erro = fabs((autoValor - autoValor_aux)/autoValor);
        autoValor_aux = autoValor;
    }
    while (erro > tolerancia);
    resultado.autoValor = autoValor;
    resultado.autoVetor = q;
    return resultado;
}

/// Função que faz a decomposição LU e encontra a inversa da matriz
/// \param matriz = Matriz original
/// \return ponteiro para a inversa da Matriz original
double** decomposicaoLU(double **matriz){
    double              **l = alocarMatriz(linhas, colunas);
    double              **u = alocarMatriz(linhas, colunas);
    double              **y = alocarMatriz(linhas, colunas);
    double **matriz_identidade = alocarMatriz(linhas, colunas);
    double  **matriz_inversa = alocarMatriz(linhas, colunas);

    //Decomposição LU -> Matriz = l * u
    //Povoamento da Matriz identidade e das Matrizes l e u
    matrizIdenditade(matriz_identidade);

    //Decomposição LU
    for(int i = 0; i < linhas; i++) {
        for(int j = 0; j < colunas; j++) {
            l[i][j] = matriz_identidade[i][j];

            if(i <= j) {
                u[i][j] = matriz[i][j];
                for(int k = 0; k < i; k++) {
                    u[i][j] -= l[i][k] * u[k][j];
                }
            }

            else {
                l[i][j] = matriz[i][j];
                for(int k = 0; k < j; k++){
                    l[i][j] -= l[i][k] * u[k][j];
                }
                l[i][j] /= u[j][j];
                u[i][j]  = 0.0;
            }
        }
    }

    //Descobrindo o y da formula para descobrir a Matriz inversa
    // (Matriz * MatrizInversa = Identidade -> l * u * MatrizInversa = Identidade -> y = u * MatrizInversa)
    for(int c = 0; c < colunas; c++) {
        for (int i = 0; i < linhas; ++i) {
            y[i][c] = matriz_identidade[i][c];
            for(int k = 0; k < i; k++) {
                y[i][c] -= l[i][k] * y[k][c];
            }
        }
    }

    //Descobrindo a MatrizInversa
    for(int c = 0; c < colunas; c++) {
        for (int i = linhas - 1; i >= 0; --i) {
            matriz_inversa[i][c] = y[i][c];
            for(int k = i + 1; k < linhas; k++) {
                matriz_inversa[i][c] -= u[i][k] * matriz_inversa[k][c];
            }
            matriz_inversa[i][c] /= u[i][i];
        }
    }
    return matriz_inversa;
}

/// Funçao que calcula o maior autovalor pelo o metodo da potencia inversa
/// \param matriz        = Matriz original
/// \param vetorInicial = Vetor chute
/// \param tolerancia     = Tolerancia para o erro
/// \return struct com o maior autovalor e o autovetor correspondente
struct vetorAutovalores potenciaInversa(double **matriz, double *vetorInicial, double tolerancia){
    struct vetorAutovalores resultado;
    struct vetorAutovalores potenciaRegular_resultado;
    double **matriz_inversa = decomposicaoLU(matriz);

    potenciaRegular_resultado   = potenciaRegular(matriz_inversa, vetorInicial, tolerancia);
    resultado.autoValor  = 1/potenciaRegular_resultado.autoValor;
    resultado.autoVetor = potenciaRegular_resultado.autoVetor;
    return resultado;
}

/// Metodo da Potencia com deslocamento
/// \param matriz
/// \param vetorInicial
/// \param tolerancia
/// \param deslocamento
/// \return uma struct com o autovalor com o deslocamento e o autovetor
struct vetorAutovalores potenciaDeslocamento(double **matriz, double *vetorInicial, double tolerancia,
                                        double deslocamento){
    double **matriz_identidade = alocarMatriz(linhas, colunas);
    double **matrizDeslocamento;
    double **matrizD_x_matrizI;

    struct vetorAutovalores resultado;
    struct vetorAutovalores potenciaInversaresultado;

    matrizIdenditade(matriz_identidade);
    matrizD_x_matrizI = matriz_x_escalar(matriz_identidade, deslocamento);
    matrizDeslocamento                 = subtrairMatriz(matriz, matrizD_x_matrizI);
    potenciaInversaresultado                   = potenciaInversa(matrizDeslocamento, vetorInicial, tolerancia);
    resultado.autoValor                  = deslocamento + potenciaInversaresultado.autoValor;
    resultado.autoVetor                 = potenciaInversaresultado.autoVetor;
    return resultado;
}

int main() {
    //Inicialização de structs
    struct vetorAutovalores           resultadopotenciaRegular;
    struct vetorAutovalores           resultadopotenciaInversa;
    struct vetorAutovalores      resultadopotenciaDeslocamento;

    double       **matriz = alocarMatriz(linhas, colunas);
    double      tolerancia = 0.00001;
    double *vetorInicial = (double*)malloc(tamanhoVetor * sizeof(double*));
    double   deslocamento = 5.0;

    double k;

    printf ("\n>>> Entrada de valores da matrícula:\n");
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

    vetorInicial[0] = 1.0;
    vetorInicial[1] = 1.0;
    vetorInicial[2] = 1.0;
    vetorInicial[3] = 1.0;
    vetorInicial[4] = 1.0;

    matriz[0][0] =  95;
    matriz[0][1] =  10;
    matriz[0][2] =  0;
    matriz[0][3] =  5;
    matriz[0][4] =  2;

    matriz[1][0] =  10;
    matriz[1][1] =  45;
    matriz[1][2] =  10;
    matriz[1][3] =  7;
    matriz[1][4] =  4;

    matriz[2][0] =  0;
    matriz[2][1] =  10;
    matriz[2][2] =  30;
    matriz[2][3] =  9;
    matriz[2][4] =  1;

    matriz[3][0] =  5;
    matriz[3][1] =  7;
    matriz[3][2] =  9;
    matriz[3][3] =  20;
    matriz[3][4] =  0.5;

    matriz[4][0] =  2;
    matriz[4][1] =  4;
    matriz[4][2] =  1;
    matriz[4][3] =  0.5;
    matriz[4][4] =  80;

    //Imprimir a matriz
    printf("\nMatriz A utilizada:\n\n");
    for (int i = 0; i < linhas ; i++) {
         for (int j = 0; j < colunas ; j++) {
            printf("%.2f\t", matriz[i][j]); 
        }
    printf("\n");
    }

    printf("\n######################## >>> Potência Regular <<< #######################\n\n");

    resultadopotenciaRegular = potenciaRegular(matriz, vetorInicial, tolerancia);
    printf("Autovalor encontrado(Potência regular): %f\n\n", resultadopotenciaRegular.autoValor);
    printf("Autovetor encontrado(Potência regular): \n\n");
    printVetor(resultadopotenciaRegular.autoVetor);

    printf("\n######################## >>> Potência Inversa <<< #######################n\n\n");

    resultadopotenciaInversa = potenciaInversa(matriz, vetorInicial, tolerancia);
    printf("Autovalor encontrado(Potência inversa): %f\n\n", resultadopotenciaInversa.autoValor);
    printf("Autovetor encontrado(Potência inversa): \n\n");
    printVetor(resultadopotenciaInversa.autoVetor);

    printf("\n######################## >>> Potência Deslocamento <<< #######################n\n\n");


    resultadopotenciaDeslocamento = potenciaDeslocamento(matriz, vetorInicial, tolerancia, deslocamento);
    printf("Autovalor encontrado(Potência deslocamento): %f\n\n", resultadopotenciaDeslocamento.autoValor);
    printf("Autovetor encontrado(Potência deslocamento): \n\n");
    printVetor(resultadopotenciaDeslocamento.autoVetor);



    k = fabs(resultadopotenciaRegular.autoValor)/fabs(resultadopotenciaInversa.autoValor);

    printf("⇥⇥⇥⇥⇥⇥⇥⇥⇥⇥⇥⇥⇥⇥⇥ >>> AP2 QUESTÃO ➀ <<< ⇤⇤⇤⇤⇤⇤⇤⇤⇤⇤⇤⇤⇤⇤⇤  : \n\nk = %f\n\n", k);




}