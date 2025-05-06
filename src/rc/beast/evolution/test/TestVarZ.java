

import rc.beast.evolution.clockmodel.MeanZCalculator;
import rc.beast.evolution.clockmodel.VarZCalculator;

public class TestVarZ {
    public static void main(String[] args) {
        double r0 = 2;
        double rt = 3.5;
        double t = 1.0;
        double psi = 0.2;
        int steps = 20000;

        double meanZ = MeanZCalculator.computeMeanZ(r0, rt, t, psi, steps);
        double varZ = VarZCalculator.computeVarZ(r0, rt, t, psi, steps);

        System.out.println("E(Z): " + meanZ);
        System.out.println("Var(Z): " + varZ);
    }


}
