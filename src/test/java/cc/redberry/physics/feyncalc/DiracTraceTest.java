package cc.redberry.physics.feyncalc;

import cc.redberry.core.TAssert;
import cc.redberry.core.indexmapping.IndexMapping;
import cc.redberry.core.indices.Indices;
import cc.redberry.core.indices.IndicesUtils;
import cc.redberry.core.parser.ParseNodeSimpleTensor;
import cc.redberry.core.parser.ParserIndices;
import cc.redberry.core.parser.preprocessor.IndicesInsertion;
import cc.redberry.core.tensor.ApplyIndexMapping;
import cc.redberry.core.tensor.Tensor;
import cc.redberry.core.tensor.Tensors;
import cc.redberry.core.transformations.expand.Expand;
import cc.redberry.core.utils.Indicator;
import cc.redberry.physics.feyncalc.context.FCC;
import junit.framework.Assert;
import org.junit.Test;

/**
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public class DiracTraceTest {
    private static Tensor parse(String str) {
        Tensor t = Tensors.parse(str, insertion);
        Indices indices = t.getIndices().getFree();
        Indices newIndices = indices.applyIndexMapping(IndexMapper.INSTANCE);
        t = ApplyIndexMapping.applyIndexMapping(t, indices.getAllIndices().copy(), newIndices.getAllIndices().copy(), new int[0]);
        return t;
    }

    private static class IndexMapper implements IndexMapping {
        public static final IndexMapper INSTANCE = new IndexMapper();
        int a = IndicesUtils.parseIndex("^a'");
        int b = IndicesUtils.parseIndex("^b'");

        @Override
        public int map(int from) {
            if (from == a)
                return b;
            return from;
        }
    }

    private static final Indicator<ParseNodeSimpleTensor> GammaIndicator = new Indicator<ParseNodeSimpleTensor>() {
        @Override
        public boolean is(ParseNodeSimpleTensor object) {
            return FCC.getGammaMatrixStringName().equals(object.name);
        }
    };

    private static final IndicesInsertion insertion = new IndicesInsertion(ParserIndices.parseSimple("^a'"), ParserIndices.parseSimple("_b'"), GammaIndicator);

    @Test
    public void test1() {
        Tensor t = parse("G_a*G_b*G_c*G_d");
        t = Expand.expand(DiracTrace.trace(t));
        Tensor expected = parse("-4*g_{ac}*g_{bd}+4*g_{ad}*g_{bc}+4*g_{ab}*g_{cd}");
        TAssert.assertEquals(t, expected);
    }

    @Test
    public void test2() {
        Tensor t = parse("G^a*G^b*G^c*G^d*G^e*G^f");
        t = Expand.expand(DiracTrace.trace(t));
        Tensor expected = parse("4*g^{af}*g^{be}*g^{cd}-4*g^{ae}*g^{bf}*g^{cd}+4*g^{ab}*g^{cd}*g^{ef}-4*g^{af}*g^{bd}*g^{ce}+4*g^{ad}*g^{bf}*g^{ce}+4*g^{ae}*g^{bd}*g^{cf}-4*g^{ad}*g^{be}*g^{cf}+4*g^{af}*g^{bc}*g^{de}-4*g^{ac}*g^{bf}*g^{de}+4*g^{ab}*g^{cf}*g^{de}-4*g^{ae}*g^{bc}*g^{df}+4*g^{ac}*g^{be}*g^{df}-4*g^{ab}*g^{ce}*g^{df}+4*g^{ad}*g^{bc}*g^{ef}-4*g^{ac}*g^{bd}*g^{ef}");
        TAssert.assertEquals(t, expected);
    }

    @Test
    public void test3() {
        Tensor[] product = {parse("G_a"), parse("G_b"), parse("G_c"), parse("G_d"), parse("G_e"), parse("G_f"), parse("G_g"), parse("G_h")};
        System.out.println(Expand.expand(DiracTrace.traceOfArray(product)));
    }

    @Test
    public void test4() {
        Tensor t = parse("G^a*G^b*G^c*G^d*G^e*G^f*G^g*G^h");
        System.out.println(Expand.expand(DiracTrace.trace(t)));
    }

    @Test
    public void test5() {
        Tensor t = parse("G^a*G^b*G^c*G^d*G^e*G^f*G^g*G^h*G^i*G^j");
        t = Expand.expand(DiracTrace.trace(t));
        Assert.assertEquals(t.size(), 945);
    }

}
