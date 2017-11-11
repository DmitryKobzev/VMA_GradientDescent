package bsu;

public class Grad_Descent {
    private final int N;
    private final double Epsilon;
    private double [][]a;
    private double [][]bMatrix;
    private double [] g;
    private double [] b;
    private double [] x;
    private double [] x_i;
    private double [] vector_nevyazki;
    private double [] nevyazki;
    private double normOfVectorNevyazki;
    Grad_Descent(){
        N=5;
        Epsilon=0.00001;
        a=new double[][]{
                {0.6444, 0.0000, -0.1683, 0.1184, 0.1973},
                {-0.0395, 0.4208, 0.0000, -0.0802, 0.0263},
                {0.0132, -0.1184, 0.7627, 0.0145, 0.0460},
                {0.0395, 0.0000, -0.0960, 0.7627, 0.0000},
                {0.0263, -0.0395, 0.1907, -0.0158, 0.5523}
        };
        bMatrix=new double[N][N];
        b=new double []{
                1.2677,
                1.6819,
                -2.3657,
                -6.5369,
                2.8351
        };
        x=new double [N];
        x_i=new double [N];
        vector_nevyazki=new double [N];
        nevyazki=new double [N];
        normOfVectorNevyazki=0;
    }
    public void printMatrix(double [][]m, int n){
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                System.out.print(m[i][j]+ " ");
            }
            System.out.println();
        }
    }
    public void printVector(double [] v , int n){
        for(int i=0;i<n;i++){
            System.out.println(v[i]+ " ");
        }
    }
    public void toSymmetricView(){
        //Транспонированная матрица
        double [][] aT=new double[N][N];
        for(int i=0;i<N;i++){
            for(int j=0;j<N;j++){
                aT[i][j]=a[j][i];
            }
        }

        System.out.println("Транспонированная матрица aT: ");
        printMatrix(aT,N);

        //умножаем транспонированную матрицу на исходную
        System.out.println("Симметричная матрица aT*a:");
        double [][] symmetricMatrix=new double[N][N];
        for(int i=0;i<N;i++){
            for(int j=0;j<N;j++){
                for(int k=0;k<N;k++){
                    symmetricMatrix[i][j] += aT[i][k]*a[k][j];
                }
            }
        }
        printMatrix(symmetricMatrix,N);

        //умножение неоднородности b на транспонированную матрицу слева
        double []f =new double [N];
       /* for(int i=0;i<N;i++){
            for(int j=0;j<N;j++){
                f[i] += aT[i][j] * b[j];
            }
        }
        for(int i=0;i<N;i++){
            b[i]=f[i];
        }*/
        System.out.println("Новая неоднородность:");
        printVector(b,N);
    }
    private double [] calc_tGradF(double [][] m,double [] x_k){
        double [] result= new double [N];
        double fraction;
        fraction=scalarMultiplication(calculateVectorNevyazki(m,x_k,N),calculateVectorNevyazki(m,x_k,N),N) / scalarMultiplication(calculateVectorNevyazki(m,x_k,N),calculateVectorNevyazki(m,x_k,N),N);
        fraction = -fraction;
        result=scalarVectorMultiplication(fraction,calculateVectorNevyazki(m,x_i,N),N);
        return result;
    }

    public void findGradientDescent(){
        int iter=0;
        while(true){
            x_i = calc_tGradF(a,x);
            for(int i=0;i<N;i++){
                x_i[i] += x[i];
            }
            iter++;
            if(isConvergesGradientDescent(x,x_i,N)== true){
                break;
            }
            for(int i=0;i<N;i++){
                x[i] = x_i[i];
            }
        }
        for(int i=0;i<N;i++){
            x[i] = x_i[i];
        }
        System.out.println("Вектор решений: ");
        printVector(x,N);
        System.out.println("Количество итераций: ");
        System.out.println(iter);
        System.out.println("Вектор невязки: ");
        nevyazki=calculateVectorNevyazki(a,x,N);
        printVector(nevyazki,N);
     //   System.out.println("Норма невязки: ");
        calculateVectorNevyazkiNorm();
    }

    private boolean isConvergesGradientDescent(double[] v1,double [] v2,int n){
        double [] v=new double[n];
        double norm=0;
        for(int i=0;i<n;i++){
            v[i]= v2[i] - v1[i];
        }
        for(int i=0;i<n;i++){
            if(Math.abs(v[i]) >= norm){
                norm=v[i];
            }
        }
        if(Math.abs(norm) < Epsilon){
            return true;
        }else{
            return false;
        }
    }

    private double scalarMultiplication(double[] v1,double [] v2,int n){
        double result = 0;
        for(int i=0;i<n;i++){
            result += v1[i] * v2[i];
        }
        return result;
    }

    private double [] matrixVectorMultiplication(double[][] m,double[] v, int n){
        double [] result=new double [n];
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                result[i] += m[i][j] * v[j];
            }
        }
        return result;
    }

    private double [] scalarVectorMultiplication(double scalar,double[] v, int n){
        double [] result=new double [n];
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                result[i] = v[i] * scalar;
            }
        }
        return result;
    }

    private double[] calculateVectorNevyazki(double [][] m,double [] x_k,int n){
        //вектор невязки
        for(int i=0;i<n;i++){
            vector_nevyazki[i]=-b[i];
            for(int j=0;j<n;j++){
                vector_nevyazki[i] += m[i][j] * x_k[j];
            }
        }
     //   System.out.println("Вектор невязки: ");
      //  printVector(vector_nevyazki,n);
        return  vector_nevyazki;
        }

    private void calculateVectorNevyazkiNorm() {
        //норма вектора невязки
        normOfVectorNevyazki = 0;
        for (int i = 0; i < N; i++) {
            if (Math.abs(vector_nevyazki[i]) >= normOfVectorNevyazki) {
                normOfVectorNevyazki = vector_nevyazki[i];
            }
        }
        System.out.println("Норма вектора невязки: ");
        System.out.println(normOfVectorNevyazki);
    }

}
