package rc.beast.evolution.clockmodel;

import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

public class GammaTransitionMatrixCalculator {

    /**
     * Container for the (alpha, beta) parameters of a Gamma distribution.
     */
    public static class AlphaBeta {
        public double alpha;
        public double beta;

        public AlphaBeta(double alpha, double beta) {
            this.alpha = alpha;
            this.beta = beta;
        }
    }

    /**
     * Computes alpha = (meanZ^2) / varZ, beta = varZ / meanZ
     * for approximating Z ~ Gamma(alpha, beta).
     */
    public static AlphaBeta computeAlphaBeta(double meanZ, double varZ) {
        if (meanZ <= 0.0 || varZ <= 0.0) {
            throw new IllegalArgumentException("meanZ and varZ must be positive.");
        }
        double alpha = (meanZ * meanZ) / varZ;
        double beta  = varZ / meanZ;
        return new AlphaBeta(alpha, beta);
    }

    /**
     * Computes the transition matrix P = (I - beta*Q)^(-alpha) via
     * eigen-decomposition and matrix power.
     *
     * @param Q      NxN rate generator matrix
     * @param alpha  gamma shape
     * @param beta   gamma scale
     * @return NxN transition probability matrix
     */
    public static double[][] computeTransitionMatrix(double[][] Q, double alpha, double beta) {
        int n = Q.length;
        RealMatrix Qmat = MatrixUtils.createRealMatrix(Q);

        // Build M = I - beta * Q
        RealMatrix I = MatrixUtils.createRealIdentityMatrix(n);
        RealMatrix M = I.subtract(Qmat.scalarMultiply(beta));

        // Eigen-decomposition of M
        EigenDecomposition ed = new EigenDecomposition(M);
        double[] realEigen = ed.getRealEigenvalues();
        double[] imagEigen = ed.getImagEigenvalues();
        RealMatrix V = ed.getV();

        // Invert V via LUDecomposition
        RealMatrix Vinv = new LUDecomposition(V).getSolver().getInverse();

        // diag(lambda^(-alpha))
        double[][] diagData = new double[n][n];
        for (int i = 0; i < n; i++) {
            if (Math.abs(imagEigen[i]) > 1e-12) {
                throw new RuntimeException(
                        "Complex eigenvalues detected. Possibly ill-conditioned Q or M."
                );
            }
            double lambda = realEigen[i];
            double powerVal = Math.pow(lambda, -alpha);
            diagData[i][i] = powerVal;
        }
        RealMatrix Dpower = MatrixUtils.createRealMatrix(diagData);

        // P = V * Dpower * V^-1
        RealMatrix Pmat = V.multiply(Dpower).multiply(Vinv);
        return Pmat.getData();
    }
}
