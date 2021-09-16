#ifndef EP2_H
#define EP2_H

#define nmax 700

// orientada por coluna
int cholcol(int n, double A[][nmax]);
int forwcol(int n, double A[][nmax], double b[]);
int backcol(int n, double A[][nmax], double b[], int trans);
// orientada por linha
int cholrow(int n, double A[][nmax]);
int forwrow(int n, double A[][nmax], double b[]);
int backrow(int n, double A[][nmax], double b[], int trans);

int lucol(int n, double A[][nmax], int p[]);
int sscol(int n, double A[][nmax], int p[], double b[]);
int lurow(int n, double A[][nmax], int p[]);
int ssrow(int n, double A[][nmax], int p[], double b[]);
int main(int argc, char *argv[]);

#endif