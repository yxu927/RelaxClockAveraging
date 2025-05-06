import rc.beast.evolution.clockmodel.MeanZCalculator;


public class TestMeanZ_SpecialCase {
    public static void main(String[] args) {
        // Example: r0=0.5, rt=2.0, t=1.5, psi=0.0
        double r0 = 0.5;
        double rt = 2.0;
        double t = 1.5;
        double psi = 0.0;
        int steps = 1000000;

        double numericMeanZ = MeanZCalculator.computeMeanZ(r0, rt, t, psi, steps);

        double v0 = Math.log(r0);
        double vt = Math.log(rt);
        // exact closed-form if psi=0
        double exactMeanZ = (Math.exp(vt) - Math.exp(v0)) / (vt - v0);

        System.out.println("Numeric MeanZ (psi=0): " + numericMeanZ);
        System.out.println("Exact MeanZ (psi=0)   : " + exactMeanZ);
    }

}