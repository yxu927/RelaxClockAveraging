open module mixture.lphybeast {
    requires beast.base;
    requires lphy.base;
    requires lphy.beast;
    requires mixture.beast;
    requires mixture.lphy;

    exports mixture.lphybeast.spi;
    exports mixture.lphybeast.tobeast;
    exports mixture.lphybeast.tobeast.generators;

    provides lphybeast.spi.LPhyBEASTMapping with mixture.lphybeast.spi.LBMixture;
}
