

import rc.beast.evolution.clockmodel.MeanZCalculator;
import rc.beast.evolution.clockmodel.VarZCalculator;

import java.util.Random;

public class TestVarZ_MonteCarlo {

    public static void main(String[] args) {
        double r0 = 2.0;
        double rt = 3.5;
        double t  = 1.0;
        double phi = 0.2; // positive variance
        int stepsForIntegration = 200; // for numeric trapezoid
        int stepsForPath = 200;       // for each MC path
        int numPaths = 30000000;        // total MC paths, increase for less noise

        // 1) Numeric integration results
        double meanZ_numerical = MeanZCalculator.computeMeanZ(r0, rt, t, phi, stepsForIntegration);
        double varZ_numerical  = VarZCalculator.computeVarZ(r0, rt, t, phi, stepsForIntegration);

        // 2) Monte Carlo simulation
        MCResults mcRes = monteCarloBridgedPaths(r0, rt, t, phi, stepsForPath, numPaths);

        // Print results
        System.out.println("=== Numeric Integration ===");
        System.out.println("E(Z): " + meanZ_numerical);
        System.out.println("Var(Z): " + varZ_numerical);
        System.out.println();
        System.out.println("=== Monte Carlo Simulation ===");
        System.out.println("E(Z): " + mcRes.meanZ);
        System.out.println("Var(Z): " + mcRes.varZ);
        System.out.println("Number of paths: " + numPaths);
    }

    /**
     * Simulate many bridged geometric Brownian paths in log-space,
     * and compute empirical mean & variance of Z = (1/t)*integral(R_s ds).
     */
    private static MCResults monteCarloBridgedPaths(
            double r0, double rt, double t, double phi,
            int stephinPath, int numPaths
    ) {
        Random rng = new Random(1234);

        double v0 = Math.log(r0);
        double vt = Math.log(rt);
        double dt = t / stephinPath;

        double sumZ = 0.0;
        double sumZ2 = 0.0;

        for (int p = 0; p < numPaths; p++) {
            double currentV = v0;
            double sumRates = 0.0; // to integrate exp(V_s)

            for (int i = 1; i <= stephinPath; i++) {
                double sPrev = (i - 1)*dt;
                double sNext = i*dt;

                // bridging drift
                double alpha = (vt - currentV)*(dt/(t - sPrev));

                // bridging variance
                double varStep = phi * dt * ((t - sNext)/(t - sPrev));
                double stdDev  = Math.sqrt(varStep);

                // random step
                double z = rng.nextGaussian();
                double nextV = currentV + alpha + stdDev*z;

                // integrate rate at midpoint
                double midV = currentV + 0.5*(nextV - currentV);
                double midRate = Math.exp(midV);

                sumRates += midRate;
                currentV = nextV;
            }

            // approximate integral from 0..t of e^{V_s} ds
            double pathIntegral = sumRates * dt;
            double pathZ = pathIntegral / t;
            sumZ += pathZ;
            sumZ2 += pathZ*pathZ;
        }

        double meanZ = sumZ / numPaths;
        double meanZ2 = sumZ2 / numPaths;
        double varZ = meanZ2 - meanZ*meanZ;

        return new MCResults(meanZ, varZ);
    }

    private static class MCResults {
        double meanZ;
        double varZ;
        MCResults(double m, double v) {
            meanZ = m;
            varZ = v;
        }
    }
}
