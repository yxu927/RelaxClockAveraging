package mixture.lphybeast.tobeast;

import beast.base.core.BEASTInterface;
import beast.base.inference.StateNode;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.spec.domain.Int;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.domain.Real;
import beast.base.spec.inference.parameter.IntScalarParam;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.RealVectorParam;
import lphy.core.model.Value;
import lphybeast.BEASTContext;

import java.util.Arrays;

public final class TypedParameterUtils {

    private TypedParameterUtils() {
    }

    public static RealVectorParam<PositiveReal> positiveRealVector(final String id,
                                                                   final double[] values,
                                                                   final boolean estimate) {
        final RealVectorParam<PositiveReal> param = new RealVectorParam<>(values, PositiveReal.INSTANCE);
        param.setID(id);
        setEstimate(param, estimate);
        param.initAndValidate();
        return param;
    }

    public static RealVectorParam<Real> realVector(final String id,
                                                   final double[] values,
                                                   final boolean estimate) {
        final RealVectorParam<Real> param = new RealVectorParam<>(values, Real.INSTANCE);
        param.setID(id);
        setEstimate(param, estimate);
        param.initAndValidate();
        return param;
    }

    public static RealScalarParam<PositiveReal> positiveRealScalar(final String id,
                                                                   final double value,
                                                                   final boolean estimate) {
        final RealScalarParam<PositiveReal> param = new RealScalarParam<>(value, PositiveReal.INSTANCE);
        param.setID(id);
        setEstimate(param, estimate);
        param.initAndValidate();
        return param;
    }

    public static RealScalarParam<Real> realScalar(final String id,
                                                   final double value,
                                                   final boolean estimate) {
        final RealScalarParam<Real> param = new RealScalarParam<>(value, Real.INSTANCE);
        param.setID(id);
        setEstimate(param, estimate);
        param.initAndValidate();
        return param;
    }

    public static IntScalarParam<Int> intScalar(final String id,
                                                final int value,
                                                final boolean estimate) {
        final IntScalarParam<Int> param = new IntScalarParam<>(value, Int.INSTANCE);
        param.setID(id);
        setEstimate(param, estimate);
        param.initAndValidate();
        return param;
    }

    public static RealVectorParam<PositiveReal> positiveRealVectorFrom(final BEASTInterface beastValue,
                                                                       final String fallbackId,
                                                                       final int targetDimension,
                                                                       final double fallbackValue,
                                                                       final boolean estimate) {
        if (beastValue instanceof RealVectorParam<?> typed && typed.getDomain() instanceof PositiveReal) {
            @SuppressWarnings("unchecked")
            final RealVectorParam<PositiveReal> positive = (RealVectorParam<PositiveReal>) typed;
            return positive;
        }

        final double[] values = valuesFromReal(beastValue, targetDimension, fallbackValue);
        final String id = idOrFallback(beastValue, fallbackId);
        return positiveRealVector(id, values, estimate);
    }

    public static RealVectorParam<Real> realVectorFrom(final BEASTInterface beastValue,
                                                       final String fallbackId,
                                                       final double[] fallbackValues,
                                                       final boolean estimate) {
        if (beastValue instanceof RealVectorParam<?> typed && typed.getDomain() instanceof Real) {
            @SuppressWarnings("unchecked")
            final RealVectorParam<Real> real = (RealVectorParam<Real>) typed;
            return real;
        }

        final double[] values = valuesFromReal(beastValue, fallbackValues.length, Double.NaN);
        final String id = idOrFallback(beastValue, fallbackId);
        return realVector(id, containsNaN(values) ? fallbackValues : values, estimate);
    }

    public static IntScalarParam<Int> intScalarFrom(final BEASTInterface beastValue,
                                                    final String fallbackId,
                                                    final int fallbackValue,
                                                    final boolean estimate) {
        if (beastValue instanceof IntScalarParam<?> typed && typed.getDomain() instanceof Int) {
            @SuppressWarnings("unchecked")
            final IntScalarParam<Int> integer = (IntScalarParam<Int>) typed;
            return integer;
        }

        final int value;
        if (beastValue instanceof IntegerParameter parameter && parameter.getDimension() >= 1) {
            value = parameter.getValue(0);
        } else if (beastValue instanceof RealParameter parameter && parameter.getDimension() >= 1) {
            value = (int) Math.round(parameter.getValue(0));
        } else {
            value = fallbackValue;
        }

        return intScalar(idOrFallback(beastValue, fallbackId), value, estimate);
    }

    public static boolean isRandom(final Value<?> value) {
        return value != null && value.isRandom();
    }

    public static double[] values(final Double[] boxed) {
        final double[] values = new double[boxed.length];
        for (int i = 0; i < boxed.length; i++) {
            values[i] = boxed[i];
        }
        return values;
    }

    public static <T extends BEASTInterface> T replaceInContext(final BEASTContext context,
                                                                final Value<?> value,
                                                                final BEASTInterface oldObject,
                                                                final T newObject) {
        if (value == null) {
            return newObject;
        }
        if (oldObject != null && oldObject != newObject) {
            context.removeBEASTObject(oldObject);
        }
        context.putBEASTObject(value, newObject);
        return newObject;
    }

    private static double[] valuesFromReal(final BEASTInterface beastValue,
                                           final int targetDimension,
                                           final double fallbackValue) {
        if (beastValue instanceof RealVectorParam<?> typed) {
            final int n = typed.size() == 1 && targetDimension > 1 ? targetDimension : typed.size();
            final double[] values = new double[n];
            for (int i = 0; i < n; i++) {
                values[i] = typed.get(typed.size() == 1 ? 0 : i);
            }
            return values;
        }
        if (beastValue instanceof RealParameter parameter) {
            final int n = parameter.getDimension() == 1 && targetDimension > 1 ? targetDimension : parameter.getDimension();
            final double[] values = new double[n];
            for (int i = 0; i < n; i++) {
                values[i] = parameter.getValue(parameter.getDimension() == 1 ? 0 : i);
            }
            return values;
        }

        final double[] values = new double[targetDimension];
        Arrays.fill(values, fallbackValue);
        return values;
    }

    private static String idOrFallback(final BEASTInterface beastValue, final String fallbackId) {
        if (beastValue != null && beastValue.getID() != null) {
            return beastValue.getID();
        }
        return fallbackId;
    }

    private static boolean containsNaN(final double[] values) {
        for (final double value : values) {
            if (Double.isNaN(value)) {
                return true;
            }
        }
        return false;
    }

    private static void setEstimate(final StateNode node, final boolean estimate) {
        node.isEstimatedInput.setValue(estimate, node);
    }
}
