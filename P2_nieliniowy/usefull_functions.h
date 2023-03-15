double* new_vector(int n) {
    double* w;
    w = (double*)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++)
        w[i] = 0.0;
    return w;
}

double** new_matrix(int W, int K) {
    double** tab;
    tab = (double**)malloc(W * sizeof(double*));
    for (int i = 0; i < W; i++)
    {
        tab[i] = (double*)malloc(K * sizeof(double));
        for (int j = 0; j < K; j++)
            tab[i][j] = 0.0;
    }
    return tab;
}

void display_vector(double* v, int n) {
    //printf("Wektor:\n");
    for (int i = 0; i < n; i++)
        printf("%lf\n", v[i]);
    printf("\n");
}

void display_matrix(double** M, int W, int K) {
    printf("Jakobian:\n");
    for (int i = 0; i < W; i++)
    {
        for (int j = 0; j < K; j++)
            printf("%.2lf\t", M[i][j]);
        printf("\n");
    }
}

void free_mat(double** M, int n) {
    for (int i = 0; i < n; ++i) {
        free(M[i]);
    }
    free(M);
}