package cc.redberry.physics.feyncalc;

import cc.redberry.core.number.Complex;
import cc.redberry.core.tensor.*;
import cc.redberry.core.transformations.ContractIndices;
import cc.redberry.core.transformations.Transformation;
import cc.redberry.core.transformations.expand.Expand;
import cc.redberry.core.transformations.expand.ExpandAll;
import cc.redberry.physics.feyncalc.context.FCC;

import java.util.ArrayList;
import java.util.List;

import static cc.redberry.core.tensor.FullContractionsStructure.*;
import static cc.redberry.core.tensor.Tensors.*;

/**
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public class DiracTrace implements Transformation {
    public static final DiracTrace DIRAC_TRACE = new DiracTrace();

    private DiracTrace() {
    }

    @Override
    public Tensor transform(Tensor t) {
        return trace(t);
    }

    public static Tensor trace(Tensor tensor) {
        tensor = ExpandAll.expandAll(tensor);
        if (!(tensor instanceof Sum))
            return traceOfProduct(tensor);
        SumBuilder sb = new SumBuilder();
        for (Tensor summand : tensor) {
            sb.put(Expand.expand(traceOfProduct(summand), ContractIndices.ContractIndices));
        }
        return sb.build();
    }

    private static Tensor traceOfProduct(Tensor tensor) {
        if (!(tensor instanceof Product))
            return multiply(Complex.FOUR, tensor);

        List<Tensor> gammaMatrices = new ArrayList<>();
        Tensor temp = tensor, current;
        for (int i = tensor.size() - 1; i >= 0; --i) {
            current = tensor.get(i);
            if (!(current instanceof SimpleTensor))
                continue;
            if (FCC.isGammaMatrix(current)) {
                gammaMatrices.add(current);
                if (temp instanceof Product)
                    temp = ((Product) temp).remove(i);
                else
                    temp = Complex.ONE;
            }
        }
        if (gammaMatrices.isEmpty())
            return multiply(Complex.FOUR, tensor);
        if (gammaMatrices.size() % 2 != 0)
            return Complex.ZERO;

        Product gammas = (Product) multiply(gammaMatrices.toArray(new Tensor[gammaMatrices.size()]));
        if (gammas.getIndices().getFree().getOfType(FCC.getDiracMatrixIndexType()).size() != 0)
            return tensor;

        return Tensors.multiply(temp, traceOfProductOfGammas(gammas));
    }

    private static Tensor traceOfProductOfGammas(Tensor tensor) {
        if (tensor instanceof SimpleTensor)
            return Complex.ZERO;
        return traceOfArray(sortProductOfMatrices((Product) tensor));
    }

    private static Tensor[] sortProductOfMatrices(Product gammas) {

        ProductContent content = gammas.getContent();
        FullContractionsStructure fs = content.getFullContractionsStructure();
        Tensor[] result = new Tensor[gammas.size()];
        result[0] = content.get(0);

        Tensor next;
        int position = 0, currentIndex = 0, nextIndex;
        long[] fromContractions;
        out:
        while (true) {
            fromContractions = fs.contractions[currentIndex];
            for (long contraction : fromContractions) {
                if (getFromIndexId(contraction) != FCC.getGammaMatricesIndicesIds().get(1))
                    continue;

                assert getToIndexId(contraction) == FCC.getGammaMatricesIndicesIds().get(0);

                next = content.get(nextIndex = getToTensorIndex(contraction));
                if (next == result[0])

                    break out;
                result[position++] = next;
                currentIndex = nextIndex;
                break;
            }
        }
        return result;
    }


    public static Tensor traceOfArray(Tensor[] product) {
        if (product.length == 2)
            return multiply(Complex.FOUR,
                    createMetricOrKronecker(product[0].getIndices().get(FCC.getLorentzIndexType(), 0),
                            product[1].getIndices().get(FCC.getLorentzIndexType(), 0)));
        if (product.length % 2 != 0)
            return Complex.ZERO;
        SumBuilder sb = new SumBuilder();
        Tensor temp;
        for (int i = 0; i < product.length - 1; ++i) {
            temp = multiply(Complex.TWO,
                    createMetricOrKronecker(product[i].getIndices().get(0), product[i + 1].getIndices().get(0)),
                    traceOfArray(subArray(product, i, i + 1)));
            if (i % 2 != 0)
                temp = negate(temp);
            sb.put(temp);
            swap(product, i, i + 1);
        }
        return multiply(Complex.ONE_HALF, sb.build());
    }

    private static Tensor[] subArray(Tensor[] array, int a, int b) {
        Tensor[] result = new Tensor[array.length - 2];
        int k = 0;
        for (int i = 0; i < array.length; ++i) {
            if (i == a || i == b)
                continue;
            result[k++] = array[i];
        }
        return result;
    }

    private static void swap(Tensor[] array, int a, int b) {
        Tensor temp = array[a];
        array[a] = array[b];
        array[b] = temp;
    }

    private static boolean isGamma(Tensor tensor, SimpleTensor gamma) {
        return tensor.getClass() == SimpleTensor.class && ((SimpleTensor) tensor).getName() == gamma.getName();
    }

}
