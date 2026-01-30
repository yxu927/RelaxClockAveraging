package mixture.lphy.spi;


import lphy.base.spi.LPhyBaseImpl;
import lphy.core.model.BasicFunction;
import lphy.core.model.GenerativeDistribution;
import mixture.lphy.evolution.auto.AutoCorrelatedClock;
import mixture.lphy.evolution.auto.AutoCorrelatedLogRates;
import mixture.lphy.evolution.auto.MixturePhyloCTMC;

import java.util.Arrays;
import java.util.List;

/**
 *
 * The provider of SPI which is an implementation of a service.
 * It requires a public no-args constructor.
 */
public class MixtureImpl extends LPhyBaseImpl {

    /**
     * Required by ServiceLoader.
     */
    public MixtureImpl() {
        //TODO print package or classes info here?
    }

    @Override
    public List<Class<? extends GenerativeDistribution>> declareDistributions() {
        return Arrays.asList(
                AutoCorrelatedLogRates.class, MixturePhyloCTMC.class

        );
    }

    @Override
    public List<Class<? extends BasicFunction>> declareFunctions() {
        return Arrays.asList(
                AutoCorrelatedClock.class


        );
    }

    public String getExtensionName() {
        return "Phylonco lphy library";
    }
}
