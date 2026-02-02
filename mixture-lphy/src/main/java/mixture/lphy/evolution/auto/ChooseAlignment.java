package mixture.lphy.evolution.auto;

import lphy.base.evolution.alignment.Alignment;
import lphy.core.model.DeterministicFunction;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorCategory;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;

import java.util.SortedMap;
import java.util.TreeMap;

public class ChooseAlignment extends DeterministicFunction<Alignment> {

    public static final String INDEX = "index";
    public static final String COMP1 = "comp1";
    public static final String COMP2 = "comp2";
    public static final String COMP3 = "comp3";

    private Value<Integer> index;
    private Value<Alignment> comp1;
    private Value<Alignment> comp2;
    private Value<Alignment> comp3; // optional

    public ChooseAlignment(
            @ParameterInfo(name = INDEX, description = "component index (0-based)") Value<Integer> index,
            @ParameterInfo(name = COMP1, description = "component 1 alignment") Value<Alignment> comp1,
            @ParameterInfo(name = COMP2, description = "component 2 alignment") Value<Alignment> comp2,
            @ParameterInfo(name = COMP3, description = "component 3 alignment", optional = true) Value<Alignment> comp3
    ) {
        this.index = index;
        this.comp1 = comp1;
        this.comp2 = comp2;
        this.comp3 = comp3;
    }

    @Override
    public SortedMap<String, Value> getParams() {
        SortedMap<String, Value> map = new TreeMap<>();
        map.put(INDEX, index);
        map.put(COMP1, comp1);
        map.put(COMP2, comp2);
        if (comp3 != null) map.put(COMP3, comp3);
        return map;
    }

    @GeneratorInfo(
            name = "ChooseAlignment",
            verbClause = "is selected from",
            narrativeName = "alignment selector",
            category = GeneratorCategory.PHYLO_LIKELIHOOD,
            description = "Deterministically selects one of the component alignments based on an index."
    )
    @Override
    public Value<Alignment> apply() {
        int k = index.value();
        Alignment chosen;

        if (k == 0) {
            chosen = comp1.value();
        } else if (k == 1) {
            chosen = comp2.value();
        } else if (k == 2 && comp3 != null) {
            chosen = comp3.value();
        } else {
            throw new IllegalArgumentException("index out of range: " + k);
        }

        if (chosen == null) {
            throw new IllegalStateException("Chosen alignment is null. Ensure component RVs are simulated before selecting.");
        }

        return new Value<>(chosen, this);
    }

    @Override
    public void setParam(String name, Value value) {
        switch (name) {
            case INDEX -> this.index = (Value<Integer>) value;
            case COMP1 -> this.comp1 = (Value<Alignment>) value;
            case COMP2 -> this.comp2 = (Value<Alignment>) value;
            case COMP3 -> this.comp3 = (Value<Alignment>) value;
            default -> throw new RuntimeException("Unknown param: " + name);
        }
    }
}
