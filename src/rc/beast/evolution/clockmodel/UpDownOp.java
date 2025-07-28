package rc.beast.evolution.clockmodel;

import beast.base.evolution.tree.Tree;
import beast.base.inference.operator.UpDownOperator;

public interface UpDownOp {

    UpDownOperator getUpDownOperator1(Tree tree);

    UpDownOperator getUpDownOperator2(Tree tree);
}
