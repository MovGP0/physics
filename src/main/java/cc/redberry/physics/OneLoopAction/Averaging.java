package cc.redberry.physics.OneLoopAction;

import cc.redberry.core.context.CC;
import cc.redberry.core.indices.IndicesBuilder;
import cc.redberry.core.indices.IndicesBuilderSorted;
import cc.redberry.core.number.ComplexElement;
import cc.redberry.core.number.NumberFraction;
import cc.redberry.core.number.RationalElement;
import cc.redberry.core.tensor.*;
import cc.redberry.transformation.Transformation;
import cc.redberry.transformation.Transformations;

public class Averaging implements Transformation {
    public static final Averaging INSTANCE = new Averaging();
    private final static SimpleTensor const_n = CC.parseSimple("n_\\mu");

    private Averaging() {
    }

    private static Tensor average(final int[] indices) {
        if (indices.length == 2)
            return CC.createMetricOrKronecker(indices[0], indices[1]);
        Sum s = new Sum();
        for (int i = 1; i < indices.length; ++i) {
            int[] suffix = new int[indices.length - 2];
            System.arraycopy(indices, 1, suffix, 0, i - 1);
            System.arraycopy(indices, i + 1, suffix, i - 1, indices.length - i - 1);
            s.add(new Product(CC.createMetricOrKronecker(indices[0], indices[i]), average(suffix)));
        }
        return s.equivalent();
    }

    @Override
    public Tensor transform(Tensor tensor) {
        if (tensor instanceof Sum) {
            TensorIterator it = tensor.iterator();
            Tensor tensorCurrent, tempResult;
            while (it.hasNext()) {
                tensorCurrent = it.next();
                tempResult = transform(tensorCurrent);
                if (TensorNumber.isZero(tempResult))
                    it.remove();
                else
                    it.set(tempResult);
            }
            return tensor.equivalent();
        }

        if (tensor instanceof Product) {
            TensorIterator it = tensor.iterator();
            int count = 0;
            Tensor current;
            IndicesBuilder ib = new IndicesBuilderSorted();
            while (it.hasNext()) {
                current = it.next();
                if (current instanceof SimpleTensor && ((SimpleTensor) current).getName() == const_n.getName()) {
                    it.remove();
                    ib.append(current);
                    ++count;
                }
            }
            if(count == 0)
                return tensor;
            if (count % 2 != 0)
                return TensorNumber.createZERO();
            count = count / 2;
            Tensor averaged = average(ib.getIndices().getAllIndices().copy());
            long factor = org.apache.commons.math.util.MathUtils.pow(2, count) * org.apache.commons.math.util.MathUtils.factorial(count + 1);
            TensorNumber number = new TensorNumber(new ComplexElement(new NumberFraction(1, factor), RationalElement.ZERO));
            averaged = Transformations.expandBrackets(averaged);
            ((Product) tensor).add(number, averaged);
        }
        return tensor.equivalent();
    }
}
