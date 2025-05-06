import java.util.Arrays;

/**
 * Demonstration code to compare:
 *   1) The "closed-form" result of Int_{0..infinity} exp(Q*z) Gamma(alpha,theta) dz
 *      i.e. (I - theta*Q)^(-alpha).
 *   2) A direct numerical integration of the same integral.
 *
 * For simplicity, everything is done for 2x2 matrices and real distinct eigenvalues.
 * This code is not production-grade; it's purely for testing the correctness of
 * the Gamma-transition matrix approach.
 */
public class GammaMatrixTest {

    public static void main(String[] args) {

        double[][] Q = {
                {-3.0, 3.0},
                { 2.0,-2.0}
        };


        double alpha = 2.0;
        double theta = 0.3;


        double[][] closedForm = computeMatrixByClosedForm(Q, alpha, theta);

        double zMax = 10.0;
        int steps = 2000;
        double[][] numeric = computeMatrixByNumericIntegration(Q, alpha, theta, zMax, steps);


        System.out.println("Closed-form result:\n" + matrixToString(closedForm));
        System.out.println("\nNumeric-integration result:\n" + matrixToString(numeric));

        double diff = compareMatrices(closedForm, numeric);
        System.out.println(String.format("\nMax element abs-difference = %.6e\n", diff));
    }


    public static double[][] computeMatrixByClosedForm(double[][] Q, double alpha, double theta) {
        // M = I - theta*Q
        double[][] M = new double[2][2];
        M[0][0] = 1.0 - theta*Q[0][0];
        M[0][1] =     - theta*Q[0][1];
        M[1][0] =     - theta*Q[1][0];
        M[1][1] = 1.0 - theta*Q[1][1];


        return matrixPower2x2(M, -alpha);
    }

    public static double[][] computeMatrixByNumericIntegration(
            double[][] Q, double alpha, double theta,
            double zMax, int steps
    ) {
        double[][] sumMat = new double[][] { {0,0},{0,0} };

        double h = zMax / steps;
        for(int i=0; i<=steps; i++){
            double z = i*h;
            double w = simpsonWeight(i, steps);

            // Evaluate e^{Q*z}
            double[][] expQz = matrixExp2x2(Q, z);

            // Evaluate gammaPDF(z)
            double pdf = gammaPDF(z, alpha, theta);

            // Accumulate
            for(int r=0;r<2;r++){
                for(int c=0;c<2;c++){
                    sumMat[r][c] += w * expQz[r][c] * pdf;
                }
            }
        }

        // multiply by step/3 if Simpson's rule
        double factor = h/3.0;
        for(int r=0;r<2;r++){
            for(int c=0;c<2;c++){
                sumMat[r][c] *= factor;
            }
        }

        return sumMat;
    }

    // Helper: gamma PDF (alpha=shape, theta=scale)
    public static double gammaPDF(double z, double alpha, double theta) {
        if(z < 0.0) return 0.0;
        // pdf(z) = 1/(Gamma(a)*theta^a) * z^(a-1) * exp(-z/theta)
        // We'll do a quick standard approach:
        double logPdf = - logGamma(alpha)
                - alpha*Math.log(theta)
                + (alpha-1.0)*Math.log(z<=0?1e-300:z)
                - (z/theta);
        return Math.exp(logPdf);
    }

    // --------------------------------------------------------------------------
    // Matrix exponent for 2x2:  e^{Q*z} using eigen-decomposition
    // We assume Q has distinct real eigenvalues.
    // --------------------------------------------------------------------------
    public static double[][] matrixExp2x2(double[][] Q, double z){

        double a = Q[0][0];
        double b = Q[0][1];
        double c = Q[1][0];
        double d = Q[1][1];

        double trace = a + d;
        double det   = a*d - b*c;
        double disc = trace*trace - 4*det;
        // assume disc>0 => distinct real eigen
        double sqrtDisc = Math.sqrt(disc);

        double lambda1 = 0.5*(trace + sqrtDisc);
        double lambda2 = 0.5*(trace - sqrtDisc);

        // V for eigen1 => (b, lambda1 - a) or something
        double[][] v1 = { {b}, {lambda1 - a} };
        double norm1 = Math.sqrt(v1[0][0]*v1[0][0] + v1[1][0]*v1[1][0]);
        if(norm1 < 1e-14){
            // try a different approach if b=0 => maybe use c
            v1[0][0] = lambda1 - d;
            v1[1][0] = c;
            norm1 = Math.sqrt(v1[0][0]*v1[0][0] + v1[1][0]*v1[1][0]);
        }
        v1[0][0]/=norm1; v1[1][0]/=norm1;

        // V for eigen2
        double[][] v2 = { {b}, {lambda2 - a} };
        double norm2 = Math.sqrt(v2[0][0]*v2[0][0] + v2[1][0]*v2[1][0]);
        if(norm2 < 1e-14){
            v2[0][0] = lambda2 - d;
            v2[1][0] = c;
            norm2 = Math.sqrt(v2[0][0]*v2[0][0] + v2[1][0]*v2[1][0]);
        }
        v2[0][0]/=norm2; v2[1][0]/=norm2;

        // Build V
        double[][] V = {
                {v1[0][0], v2[0][0]},
                {v1[1][0], v2[1][0]}
        };
        double[][] Vinv = invert2x2(V);

        // diag( e^(lambda1*z), e^(lambda2*z) )
        double e1 = Math.exp(lambda1*z);
        double e2 = Math.exp(lambda2*z);
        double[][] diagE = {
                { e1, 0 },
                { 0,  e2}
        };

        // e^{Qz} = V * diagE * V^-1
        double[][] tmp = mul2x2(V, diagE);
        double[][] out = mul2x2(tmp, Vinv);
        return out;
    }


    public static double[][] matrixPower2x2(double[][] M, double p) {
        // 1) find eigen, M= V diag(l1,l2) V^-1
        // 2) M^p= V diag(l1^p, l2^p) V^-1
        double a = M[0][0];
        double b = M[0][1];
        double c = M[1][0];
        double d = M[1][1];

        double trace = a + d;
        double det   = a*d - b*c;
        double disc = trace*trace - 4.0*det;
        if(disc<0) {
            throw new RuntimeException("matrixPower2x2: negative discriminant => complex eigen, not handled!");
        }
        double sqrtDisc = Math.sqrt(disc);
        double lambda1 = 0.5*(trace + sqrtDisc);
        double lambda2 = 0.5*(trace - sqrtDisc);

        // eigenvector 1
        double[][] v1 = { {b}, {lambda1 - a} };
        double norm1 = Math.sqrt(v1[0][0]*v1[0][0] + v1[1][0]*v1[1][0]);
        if(norm1<1e-14){
            v1[0][0] = lambda1 - d;
            v1[1][0] = c;
            norm1 = Math.sqrt(v1[0][0]*v1[0][0] + v1[1][0]*v1[1][0]);
        }
        v1[0][0]/=norm1; v1[1][0]/=norm1;

        // eigenvector 2
        double[][] v2 = { {b}, {lambda2 - a} };
        double norm2 = Math.sqrt(v2[0][0]*v2[0][0] + v2[1][0]*v2[1][0]);
        if(norm2<1e-14){
            v2[0][0] = lambda2 - d;
            v2[1][0] = c;
            norm2 = Math.sqrt(v2[0][0]*v2[0][0] + v2[1][0]*v2[1][0]);
        }
        v2[0][0]/=norm2; v2[1][0]/=norm2;

        double[][] V = {
                {v1[0][0], v2[0][0]},
                {v1[1][0], v2[1][0]}
        };
        double[][] Vinv = invert2x2(V);

        // diag(l1^p, l2^p)
        if(lambda1<=0 || lambda2<=0) {
            // might cause complex or infinite if p is negative
            // or zero eigen => not invertible, etc.
            // you might want extra checks here
        }
        double l1p = Math.pow(lambda1, p);
        double l2p = Math.pow(lambda2, p);

        double[][] diagPow = { {l1p, 0},{0, l2p} };

        double[][] tmp = mul2x2(V, diagPow);
        double[][] out = mul2x2(tmp, Vinv);
        return out;
    }


    public static double[][] invert2x2(double[][] A){
        double a = A[0][0], b = A[0][1], c = A[1][0], d = A[1][1];
        double det = a*d - b*c;
        if(Math.abs(det)<1e-14){
            throw new RuntimeException("invert2x2: det~0 => not invertible");
        }
        double invdet = 1.0/det;
        return new double[][] {
                { d*invdet, -b*invdet},
                {-c*invdet,  a*invdet}
        };
    }

    public static double[][] mul2x2(double[][] A, double[][] B){
        double[][] C = new double[2][2];
        C[0][0] = A[0][0]*B[0][0] + A[0][1]*B[1][0];
        C[0][1] = A[0][0]*B[0][1] + A[0][1]*B[1][1];
        C[1][0] = A[1][0]*B[0][0] + A[1][1]*B[1][0];
        C[1][1] = A[1][0]*B[0][1] + A[1][1]*B[1][1];
        return C;
    }


    private static double simpsonWeight(int i, int n) {
        // standard 1D Simpson on [0..n]
        // w(0)=w(n)=1, odd =>4, even =>2
        if(i==0 || i==n) return 1.0;
        if((i % 2)==1) return 4.0;
        return 2.0;
    }

    public static double compareMatrices(double[][] A, double[][] B){
        double maxDiff=0.0;
        for(int r=0;r<2;r++){
            for(int c=0;c<2;c++){
                double diff = Math.abs(A[r][c] - B[r][c]);
                if(diff>maxDiff) maxDiff=diff;
            }
        }
        return maxDiff;
    }


    public static String matrixToString(double[][] M){
        return String.format("[[%.6f, %.6f],\n [%.6f, %.6f]]",
                M[0][0], M[0][1], M[1][0], M[1][1]);
    }


    public static double logGamma(double x){
        // Quick Lanczos approximation, for x>0
        // constants
        double[] cof = {
                76.18009172947146,
                -86.50532032941677,
                24.01409824083091,
                -1.231739572450155,
                0.1208650973866179e-2,
                -0.5395239384953e-5
        };
        double xx = x;
        double tmp = x + 5.5;
        tmp -= (xx + 0.5)*Math.log(tmp);
        double ser = 1.000000000190015;
        for (int j=0; j<6; j++){
            xx += 1.0;
            ser += cof[j]/xx;
        }
        return -tmp + Math.log(2.5066282746310005*ser/x);
    }
}
