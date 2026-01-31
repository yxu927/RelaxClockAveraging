
module mixture.lphy {
    requires transitive lphy.base;
    requires jdk.jfr;


    exports mixture.lphy.evolution.auto;

    // declare what service interface the provider intends to use
    uses lphy.core.spi.Extension;
    provides lphy.core.spi.Extension with mixture.lphy.spi.MixtureImpl;
}