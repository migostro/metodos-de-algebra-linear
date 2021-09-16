#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "ep2.h"

#define PARTE 1 // define qual parte vai fazer os testes

const double eps = 0.00000001;

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        printf("Está faltando um parâmetro!\nExecute como:\n./ep2 <nome_arquivo.dat>\n");
    }
    
    // B e b_linha são apenas matrizes auxiliares para reiniciar a valores iniciais a matriz A e b para executar o proximo teste
    double A[nmax][nmax], B[nmax][nmax], b[nmax], b_linha[nmax];
    double aux, tempo;
    int n, p[nmax];
    FILE * arquivo;
    struct timeval comeco, fim;

    arquivo = fopen(argv[1],"r");

    fscanf(arquivo, "%d", &n);
    
    int j, k;
    for (int i = 0; i < n; i++)
    {
        for (int ii = 0; ii < n; ii++)
        {
            fscanf(arquivo, "%d %d %lf", &j, &k, &aux);
            A[k][j] = aux;
            B[k][j] = aux;
        }
    }

    for (int i = 0; i < n; i++)
    {
        fscanf(arquivo, "%d %lf", &k, &aux);
        b[k] = aux;
        b_linha[k] = aux;
    }

    /* TESTES PARTE 1 */
    if (PARTE == 1)
    {  
        gettimeofday(&comeco, NULL);
        cholcol(n, A);
        gettimeofday(&fim, NULL);
        tempo = (double) fim.tv_sec - comeco.tv_sec + (double) (fim.tv_usec - comeco.tv_usec)/1000000;
        printf("cholcol: %lf\n", tempo);

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                A[i][j] = B[i][j];
            }
            b[i] = b_linha[i]; 
        }

        gettimeofday(&comeco, NULL);
        forwcol(n, A, b);
        gettimeofday(&fim, NULL);
        tempo = (double) fim.tv_sec - comeco.tv_sec + (double) (fim.tv_usec - comeco.tv_usec)/1000000;
        printf("forwcol: %lf\n", tempo);

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                A[i][j] = B[i][j];
            }
            b[i] = b_linha[i]; 
        }

        gettimeofday(&comeco, NULL);
        backcol(n, A, b, 0);
        gettimeofday(&fim, NULL);
        tempo = (double) fim.tv_sec - comeco.tv_sec + (double) (fim.tv_usec - comeco.tv_usec)/1000000;
        printf("backcol: %lf\n", tempo);

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                A[i][j] = B[i][j];
            }
            b[i] = b_linha[i]; 
        }

        gettimeofday(&comeco, NULL);
        cholrow(n, A);
        gettimeofday(&fim, NULL);
        tempo = (double) fim.tv_sec - comeco.tv_sec + (double) (fim.tv_usec - comeco.tv_usec)/1000000;
        printf("cholrow: %lf\n", tempo);

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                A[i][j] = B[i][j];
            }
            b[i] = b_linha[i]; 
        }

        gettimeofday(&comeco, NULL);
        forwrow(n, A, b);
        gettimeofday(&fim, NULL);
        tempo = (double) fim.tv_sec - comeco.tv_sec + (double) (fim.tv_usec - comeco.tv_usec)/1000000;
        printf("forwrow: %lf\n", tempo);

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                A[i][j] = B[i][j];
            }
            b[i] = b_linha[i]; 
        }

        gettimeofday(&comeco, NULL);
        backrow(n, A, b, 0);
        gettimeofday(&fim, NULL);
        tempo = (double) fim.tv_sec - comeco.tv_sec + (double) (fim.tv_usec - comeco.tv_usec)/1000000;
        printf("backrow: %lf\n", tempo);
    }
    else /* TESTES PARTE 2 */
    {
        gettimeofday(&comeco, NULL);
        lucol(n, A, p);
        gettimeofday(&fim, NULL);
        tempo = (double) fim.tv_sec - comeco.tv_sec + (double) (fim.tv_usec - comeco.tv_usec)/1000000;
        printf("lucol: %lf\n", tempo);

        gettimeofday(&comeco, NULL);
        sscol(n, A, p, b);
        gettimeofday(&fim, NULL);
        tempo = (double) fim.tv_sec - comeco.tv_sec + (double) (fim.tv_usec - comeco.tv_usec)/1000000;
        printf("sscol: %lf\n", tempo);

        // restaura a matriz e o vetor ao estado inicial
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                A[i][j] = B[i][j];
            }
            b[i] = b_linha[i]; 
        }
        
        gettimeofday(&comeco, NULL);
        lurow(n, A, p);
        gettimeofday(&fim, NULL);
        tempo = (double) fim.tv_sec - comeco.tv_sec + (double) (fim.tv_usec - comeco.tv_usec)/1000000;
        printf("lurow: %lf\n", tempo);

        gettimeofday(&comeco, NULL);
        ssrow(n, A, p, b);
        gettimeofday(&fim, NULL);
        tempo = (double) fim.tv_sec - comeco.tv_sec + (double) (fim.tv_usec - comeco.tv_usec)/1000000;
        printf("ssrow: %lf\n", tempo);
    }
    
    return 0;
}

/******* ORIENTADA POR COLUNA *******/

int cholcol(int n, double A[][nmax])
{
    for (int i = 0; i < n; i++)
    {
        if (A[i][i] == 0)
        {
            return -1;
        }

        // calcula o R_ii
        for (int k = 0; k < i; k++)
        {
            A[i][i] = A[i][i] - pow(A[k][i], 2);
        }
        A[i][i] = sqrt(A[i][i]);

        // calcula para o restante dos R_ij para j > i
        for (int j = i+1; j < n; j++)
        {
            for (int k = 0; k < i-1; k++)
            {
                A[i][j] = A[i][j] - A[k][i]*A[k][j];
            }
            A[i][j] /= A[i][i];
        }
    }

    return 0;
}

int forwcol(int n, double A[][nmax], double b[])
{
    for (int j = 0; j < n; j++)
    {
        // verifica se a matriz é singular
        if (A[j][j] == 0)
        {
            return -1;
        }

        for (int i = 0; i < j; i++)
        {
            b[i] = b[i] - A[i][j]*b[j];
        }

        b[j] = b[j]/A[j][j];
    }

    return 0;
}

int backcol(int n, double A[][nmax], double b[], int trans)
{
    for (int j = n-1; j >= 0; j--)
    {
        // verifica se a matriz é singular
        if (A[j][j] == 0)
        {
            return -1;
        }

        for (int i = n-1; i > j; i--)
        {
            b[i] = b[i] - A[i][j]*b[j];
        }
        
        b[j] = b[j]/A[j][j];
    }

    return 0;
}

/******* ORIENTADA POR LINHA *******/

int cholrow(int n, double A[][nmax])
{
    for (int i = 0; i < n; i++)
    {
        if (A[i][i] == 0)
        {
            return -1;
        }

        // calcula o R_ii
        for (int k = 0; k < i; k++)
        {
            A[i][i] = A[i][i] - pow(A[k][i], 2);
        }
        A[i][i] = sqrt(A[i][i]);

        // calcula para o restante dos R_ij para j > i
        for (int k = 0; k < i-1; k++)
        {
            for (int j = i+1; j < n; j++)
            {
                A[i][j] = (A[i][j] - A[k][i]*A[k][j]) / A[i][i];
            }
        }
    }

    return 0;
}

int forwrow(int n, double A[][nmax], double b[])
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < i; j++)
        {
            b[i] = b[i] - A[i][j]*b[j];
        }
        // verifica se a matriz é singular
        if (A[i][i] == 0)
        {
            return -1;
        }
        b[i] = b[i]/A[i][i];
    }

    return 0;
}

int backrow(int n, double A[][nmax], double b[], int trans)
{
    for (int i = n-1; i >= 0; i--)
    {
        for (int j = n-1; j > i; j--)
        {
            b[i] = b[i] - A[i][j]*b[j];
        }
        // verifica se a matriz é singular
        if (A[i][i] == 0)
        {
            return -1;
        }
        b[i] = b[i]/A[i][i];
    }
    return 0;
}

int lucol(int n, double A[][nmax], int p[]){
    double max = -1, aux;
    int m;
    for (int i = 0; i < n-1; i++)
    {
        // acha o maior valor da coluna
        for (int j = i; j < n; j++)
        {
            if (abs(A[j][i]) > max)
            {
                max = abs(A[j][i]);
                m = j;
            }
            
        }
        // verifica se A é singular
        if (max != 0)
        {
            p[i] = m;
            
            // verifica se tem a necessidade de trocar as linhas
            if (m != i)
            {
                // faz a troca de linha de matriz A
                for (int j = i; j < n; j++)
                {
                    aux = A[i][j];
                    A[i][j] = A[i][m];
                    A[i][m] = aux;
                }
                
                // faz a troca de posição do vetor p
                aux = p[i];
                p[i] = p[m];
                p[m] = aux;
            }

            // acha as os multiplos necessarios para zerar ao subtrair as linhas
            for (int j = i+1; j < n; j++)
            {
                A[j][i] = A[j][i]/A[i][i];
            }
            
            // subtrai a k vezes a linha i da linha j
            // sendo k o numero de vezes que necessario para zerar a linha
            // calcula em operações por colunas
            for (int j = i+1; j < n; j++)
            {
                for (int k = i+1; k < n; k++)
                {
                    A[k][j] = A[k][j] - A[k][i]*A[i][j];
                }
            }
        }
        else
        {
            return -1;
        }
        
    }
    // verifica se A é singular
    if (A[n-1][n-1] == 0)
    {
        return -1;
    }
    else
    {
        p[n-1] = n-1;
    }
    
    return 0;
}

int sscol(int n, double A[][nmax], int p[], double b[]){
    int m;
    double aux;
    // faz as trocas de itens correspondentes as trocas de linha da matriz A
    for (int i = 0; i < n-1; i++)
    {
        m = p[i];

        aux = b[i];
        b[i] = b[m];
        b[m] = aux;
    }
    
    // primeira parte do código que é orientada a coluna
    for (int j = 0; j < n-1; j++)
    {
        for (int i = j+1; i < n; i++)
        {
            b[i] = b[i] - A[i][j]*b[j];
        }
    }
    
    // segunda parte do código que é orientada a coluna
    for (int j = n-1; j >= 0; j++)
    {
        // verifica se a matriz é singular
        if (A[j][j] == 0)
        {
            return -1;
        }// calcula a solução do problema
        b[j] = b[j]/A[j][j];
        for (int i = 0; i < j; i++)
        {
            b[i] = b[i] - A[i][j]*b[j];
        }
    }
    
    return 0;
}

int lurow(int n, double A[][nmax], int p[]){
    double max = -1, aux;
    int m;
    for (int i = 0; i < n-1; i++)
    {
        // acha o maior valor da coluna
        for (int j = i; j < n; j++)
        {
            if (abs(A[j][i]) > max)
            {
                max = abs(A[j][i]);
                m = j;
            }
            
        }
        // verifica se A é singular
        if (max != 0)
        {
            p[i] = m;
            
            // verifica se tem a necessidade de trocar as linhas
            if (m != i)
            {
                // faz a troca de linha da matriz A
                for (int j = 0; j < n; j++)
                {
                    aux = A[i][j];
                    A[i][j] = A[i][m];
                    A[i][m] = aux;
                }
                // faz a troca de posição do vetor p
                aux = p[i];
                p[i] = p[m];
                p[m] = aux;
            }

            // acha as os multiplos necessarios para zerar ao subtrair as linhas
            for (int j = i+1; j < n; j++)
            {
                A[j][i] = A[j][i]/A[i][i];
            }
            
            // subtrai a k vezes a linha i da linha j
            // sendo k o numero de vezes que necessario para zerar a linha
            // calcula em operações por linhas
            for (int k = i+1; k < n; k++)
            {
                for (int j = i+1; j < n; j++)
                {
                    A[k][j] = A[k][j] - A[k][i]*A[i][j];
                }
            }
        }
        else
        {
            return -1;
        }
        
    }
    // verifica se A é singular
    if (A[n-1][n-1] == 0)
    {
        return -1;
    }
    else
    {
        p[n-1] = n-1;
    }
    
    return 0;
}

int ssrow(int n, double A[][nmax], int p[], double b[]){
    int m;
    double aux;
    // faz as trocas de itens correspondentes as trocas de linha da matriz A
    for (int i = 0; i < n-1; i++)
    {
        m = p[i];

        aux = b[i];
        b[i] = b[m];
        b[m] = aux;
    }
    
    // calcula a solução Pb orientado a linha
    for (int j = 0; j < n-1; j++)
    {
        for (int i = j+1; i < n; i++)
        {
            b[i] = b[i] - A[i][j]*b[j];
        }
    }
    
    // segunda parte do código que é orientada a linha
    for (int i = n-1; i >= 0; i++)
    {
        // verifica se a matriz é singular
        if (A[i][i] == 0)
        {
            return -1;
        }
        // calcula a solução do problema
        b[i] = b[i]/A[i][i];
        for (int j = 0; j < i; j++)
        {
            b[i] = b[i] - A[i][j]*b[j];
        }
    }

    return 0;
}