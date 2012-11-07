package cc.redberry.physics.feyncalc;

import cc.redberry.core.TAssert;
import cc.redberry.core.indices.IndexType;
import cc.redberry.core.parser.ParseNodeSimpleTensor;
import cc.redberry.core.parser.ParserIndices;
import cc.redberry.core.parser.preprocessor.IndicesInsertion;
import cc.redberry.core.tensor.Product;
import cc.redberry.core.tensor.Tensor;
import cc.redberry.core.tensor.Tensors;
import cc.redberry.core.utils.Indicator;
import org.junit.Test;

import static cc.redberry.physics.feyncalc.InverseProductOfMatrices.inverseProductsOfMatrices;

/**
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public class InverseProductOfMatricesTest {
    private static final Indicator<ParseNodeSimpleTensor> matricesIndicator
            = new Indicator<ParseNodeSimpleTensor>() {
        @Override
        public boolean is(ParseNodeSimpleTensor t) {
            if (t.name.length() != 1)
                return false;
            return Character.isUpperCase(t.name.charAt(0));
        }
    };

    private static final IndicesInsertion indicesInsertion = new IndicesInsertion(
            ParserIndices.parseSimple("^a'"),
            ParserIndices.parseSimple("_b'"), matricesIndicator);

    private static Tensor parseMatrices(String str) {
        return Tensors.parse(str, indicesInsertion);
    }

    @Test
    public void test1() {
        Product p = (Product) parseMatrices("A*B*C");
        TAssert.assertEquals(inverseProductsOfMatrices(p, IndexType.LatinLower1),
                parseMatrices("C*B*A"));
        TAssert.assertNotEquals(inverseProductsOfMatrices(p, IndexType.LatinLower1),
                parseMatrices("C*A*B"));
    }


    @Test
    public void test2() {
        Product p = (Product) parseMatrices("a*A*B*C");
        TAssert.assertEquals(inverseProductsOfMatrices(p, IndexType.LatinLower1),
                parseMatrices("C*a*B*A"));
        TAssert.assertNotEquals(inverseProductsOfMatrices(p, IndexType.LatinLower1),
                parseMatrices("C*A*a*B"));
    }


    @Test
    public void test3() {
        Product p = (Product) parseMatrices("a*A_i*B^ip*C_p");
        TAssert.assertEquals(inverseProductsOfMatrices(p, IndexType.LatinLower1),
                parseMatrices("C_p*a*B^ip*A_i"));
        TAssert.assertNotEquals(inverseProductsOfMatrices(p, IndexType.LatinLower1),
                parseMatrices("C_p*A_i*a*B^ip"));
    }

    @Test
    public void test4() {
        Tensor t = parseMatrices("a*A_i*B^ip*C_p - Sin[x]*D*E");
        TAssert.assertEquals(inverseProductsOfMatrices(t, IndexType.LatinLower1),
                parseMatrices("C_p*a*B^ip*A_i - Sin[x]*E*D"));
        TAssert.assertNotEquals(inverseProductsOfMatrices(t, IndexType.LatinLower1),
                parseMatrices("C_p*A_i*a*B^ip - Sin[x]*D*E"));
        TAssert.assertNotEquals(inverseProductsOfMatrices(t, IndexType.LatinLower1),
                parseMatrices("a*A_i*B^ip*C_p - Sin[x]*D*E"));
    }

}
