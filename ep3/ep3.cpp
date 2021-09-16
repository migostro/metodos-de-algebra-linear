#include <iostream>
#include <cmath>
#include <sys/time.h>

using namespace std;

// precisão da máquina = 10^-8
#define epsilon 0.0000000001

// Faz a decomposição QR por refletores com pivoteamento de colunas
// Resolve tanto sistemas lineares com posto completo quanto incompleto
void QR_refletores (int n, int m, double ** A)
{

    double Beta = -1, gama, r, maior_norma = -1, aux = -1, norma_A, precisao, norma_x, norma_k, sum;
    double x_aux[n], x[n], v[n];
    int maior_indice;

/***************************************
    // encontra o maior elemento de A para, adiante, usar no scalling
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            if (A[i][j] > Beta)
            {
                Beta = A[i][j];
            }
        }
    }
    cout << Beta << endl;

    // aplica o scalling para evitar overflows e underflows
    // E aproveita para encontrar a Norma de A para considerar oq será ou não zero
    // calculando a norma antes mesmo de aplicar o scalling para achar a norma verdadeira e não a dividida por max{a_ij}
    for (int i = 0; i < n; i++)
    {
        x_aux[i] = 0;
        for (int j = 0; j < m; j++)
        {
            x_aux[i] += A[i][j]*x[i]; // norma de A
            A[i][j] /= Beta; // scalling
        }

        if (x_aux[i] > norma_A)
            norma_A = x_aux[i];

        if (x[i] > maior_norma)
            maior_norma = x[i];
    }

    for (int k = 0; k < m; k++)
    {
        maior_norma = -1;
        maior_indice = k;

        // procura a coluna de maior norma para fazer o pivoteamento das colunas
        // orientado a coluna já que aqui não há como fazer diferente
        for (int i = k; i < m; i++)
        {
            //aux_norma = 0;
            for (int j = 0; j < n; j++)
            {
                if (abs(A[j][i]) > maior_norma)
                {
                    maior_norma = abs(A[j][i]);
                    maior_indice = i;
                }
            }
        }

        // troca as colunas caso maior indice diferente de l
        if (maior_indice != k)
        {
            for (int i = k; i < n; i++)
            {
                aux = A[i][k];
                A[i][k] = A[i][maior_indice];
                A[i][maior_indice] = aux;
            }
        }

        norma_x = 0;
        for (int i = k; i < n; i++)
        {
            v[k] = A[i][k];

            if(A[i][k] > norma_x)
            {
                norma_x += pow(A[i][k], 2);
            }

        }
        norma_x = sqrt(norma_x);

        if (v[0] > 0)
        {
            v[k] = v[k] + norma_x;
        }
        else
        {
            v[k] = v[k] - norma_x;
        }
        //v[k] = v[k] + (A[0][0]/abs(A[0][0]))*norma_x;


        norma_k = 0;
        for (int i = k; i < n; i++)
        {
            norma_k += pow(v[k], 2);
        }

        norma_k = sqrt(norma_k);
        
        for (int i = k; i < n; i++)
        {
            v[k] = v[k]/norma_k;
        }

        // Calcula a multiplicação vA 
        for (int i = k; i < n; i++)
        {
            for (int j = k; j < m; j++)
            {
                x_aux[i] = v[i]*A[i][j];
            }
        }

        for (int i = k; i < n; i++)
        {
            sum = x_aux[i]*v[i];
        }
        
        
        for (int i = k; i < n; i++)
        {
            for (int j = k; j < m; j++)
            {
                A[i][j] = A[i][j] - 2*sum;
            }
        }
    }
    
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            A[i][j] = A[i][j]*Beta;
        }
    }
    

******************************************************************************/
    // encontra o maior elemento de A para, adiante, usar no scalling
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            if (A[i][j] > Beta)
            {
                Beta = A[i][j];
            }
        }
    }
    cout << Beta << endl;

    // aplica o scalling para evitar overflows e underflows
    // E aproveita para encontrar a Norma de A para considerar oq será ou não zero
    // calculando a norma antes mesmo de aplicar o scalling para achar a norma verdadeira e não a dividida por max{a_ij}
    for (int i = 0; i < n; i++)
    {
        x_aux[i] = 0;
        for (int j = 0; j < m; j++)
        {
            x_aux[i] += A[i][j]*x[i]; // norma de A
            A[i][j] /= Beta; // scalling
        }

        if (x_aux[i] > norma_A)
            norma_A = x_aux[i];

        if (x[i] > maior_norma)
            maior_norma = x[i];
    }
    
    // resolve ||Ax||/||x||, encontrando a norma de A
    norma_A /= maior_norma;

    precisao = norma_A*epsilon;

    for (int k = 0; k < n-1; k++)
    {
        maior_norma = -1;
        maior_indice = k;

        // procura a coluna de maior norma para fazer o pivoteamento das colunas
        // orientado a coluna já que aqui não há como fazer diferente
        for (int i = k; i < m; i++)
        {
            //aux_norma = 0;
            for (int j = 0; j < n; j++)
            {
                if (abs(A[j][i]) > maior_norma)
                {
                    maior_norma = abs(A[j][i]);
                    maior_indice = i;
                }
            }
        }

        // troca as colunas caso maior indice diferente de l
        if (maior_indice != k)
        {
            for (int i = k; i < n; i++)
            {
                aux = A[i][k];
                A[i][k] = A[i][maior_indice];
                A[i][maior_indice] = aux;
            }
        }

        // verifica se o maior elemento da diagonal pricipal é numericamente 0
        // caso positivo acabamos o métodos pois as linhas de k até n são todos zeros
        if (abs(A[k][k]) < epsilon)
        {
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < m; j++)
                {
                    A[i][j] = A[i][j]*Beta;
                }
            }
            /*cout << "A " << A[k][k] << endl;
            cout << "precisão " << precisao << endl;
            cout << "SAIU " << k << endl;*/
            return;
        }

        // encontra gama, r e u tal que Q = (I-gama*u*u^t), sendo t o transposto do vetor u
        // e u sendo armazenado abaixo de do vetor x

        r = 0;



        // calcula r da coluna de A
        for (int i = k; i < n; i++)
        {
            r += pow(A[i][k], 2);
        }

        r = sqrt(r);
        
        if (A[k][k] < 0)
            r = -r;
        
        A[k][k] = r + A[k][k];
        
        gama = A[k][k]/r;

        // armazena u embaixo de x
        for (int i = k+1; i < n; i++)
            A[i][k] = A[i][k]/A[k][k];
        
        A[k][k] = 1;

        
        // calcula a matriz (n-k)x(n-k) calculando I - gama_k*u_k*u_k^t sendo t transposto
        //v = gama*u;
        for (int i = k; i < n; i++)
        {
            v[i] = gama*A[i][k];
        }
        //v = v*B;
        for (int i = k; i < n; i++)
        {
            for (int j = k; j < m; j++)
            {
                v[i] = v[i]*A[i][j];
            }

            aux = v[i]*A[i][k];
        }
        //B = B - uv;
        for (int i = k; i < n; i++)
        {
            for (int j = k+1; j < m; j++)
            {
                A[i][j] = A[i][j] - aux;
            }
        }
        
        //A[k][k] = -r;
    }

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            A[i][j] = A[i][j]*Beta;
        }
    }
//****************************************** 
}

int main(int argc, char *argv[])
{
    int n, m;
    cout << "Coloque o tamanho de n e m:" << endl;
    cin >> n >> m;

    double ** A;

    struct timeval tempo_inicio, tempo_fim;

    A = (double **) malloc(n*sizeof(double*));
    for (int i = 0; i < n; i++)
    {
        A[i] = (double *) malloc(m*sizeof(double));
    }

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            A[i][j] = rand();
        }
    }

    /*for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            cout << A[i][j] << " ";
        }
        cout << endl;
    }*/

    gettimeofday(&tempo_inicio, NULL);
    QR_refletores(n,m,A);
    gettimeofday(&tempo_fim, NULL);

    cout << endl << "TEMPO: " << (double) (tempo_fim.tv_usec - tempo_inicio.tv_usec) / 1000000 + (double) (tempo_fim.tv_sec - tempo_inicio.tv_sec) << endl;
    /*for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            cout << A[i][j] << " ";
        }
        cout << endl;
    }*/
    
}