package cc.redberry.physics.feyncalc;

import cc.redberry.core.TAssert;
import cc.redberry.core.context.CC;
import cc.redberry.core.indices.IndexType;
import cc.redberry.core.parser.preprocessor.GeneralIndicesInsertion;
import cc.redberry.core.tensor.Tensor;
import cc.redberry.core.transformations.expand.Expand;
import junit.framework.Assert;
import org.junit.Test;

import static cc.redberry.core.tensor.Tensors.parse;
import static cc.redberry.core.tensor.Tensors.parseSimple;

/**
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public class DiracTraceTest {
    @Test
    public void test1() {
        GeneralIndicesInsertion indicesInsertion = new GeneralIndicesInsertion();
        CC.current().getParseManager().defaultParserPreprocessors.add(indicesInsertion);
        indicesInsertion.addInsertionRule(parseSimple("G^a'_b'a"), IndexType.LatinLower1);

        Tensor t = parse("Tr[G_a*G_b*G_c*G_d]");
        t = Expand.expand(DiracTrace.trace(t));
        Tensor expected = parse("-4*g_{ac}*g_{bd}+4*g_{ad}*g_{bc}+4*g_{ab}*g_{cd}");
        TAssert.assertEquals(t, expected);
    }

    @Test
    public void test2() {
        GeneralIndicesInsertion indicesInsertion = new GeneralIndicesInsertion();
        CC.current().getParseManager().defaultParserPreprocessors.add(indicesInsertion);
        indicesInsertion.addInsertionRule(parseSimple("G^a'_b'a"), IndexType.LatinLower1);

        Tensor t = parse("Tr[G^a*G^b*G^c*G^d*G^e*G^f]");
        t = Expand.expand(DiracTrace.trace(t));
        Tensor expected = parse("4*g^{af}*g^{be}*g^{cd}-4*g^{ae}*g^{bf}*g^{cd}+4*g^{ab}*g^{cd}*g^{ef}-4*g^{af}*g^{bd}*g^{ce}+4*g^{ad}*g^{bf}*g^{ce}+4*g^{ae}*g^{bd}*g^{cf}-4*g^{ad}*g^{be}*g^{cf}+4*g^{af}*g^{bc}*g^{de}-4*g^{ac}*g^{bf}*g^{de}+4*g^{ab}*g^{cf}*g^{de}-4*g^{ae}*g^{bc}*g^{df}+4*g^{ac}*g^{be}*g^{df}-4*g^{ab}*g^{ce}*g^{df}+4*g^{ad}*g^{bc}*g^{ef}-4*g^{ac}*g^{bd}*g^{ef}");
        TAssert.assertEquals(t, expected);
    }

    @Test
    public void test3() {
        Tensor[] product = {parse("G_a"), parse("G_b"), parse("G_c"), parse("G_d"), parse("G_e"), parse("G_f"), parse("G_g"), parse("G_h")};
        System.out.println(Expand.expand(DiracTrace.traceOfArray(product, IndexType.LatinLower)));
    }

    @Test
    public void test4() {
        GeneralIndicesInsertion indicesInsertion = new GeneralIndicesInsertion();
        CC.current().getParseManager().defaultParserPreprocessors.add(indicesInsertion);
        indicesInsertion.addInsertionRule(parseSimple("G^a'_b'a"), IndexType.LatinLower1);

        Tensor t = parse("Tr[G^a*G^b*G^c*G^d*G^e*G^f*G^g*G^h]");
        t = Expand.expand(DiracTrace.trace(t));
        TAssert.assertEquals(t, "-4*g^{eh}*g^{ac}*g^{fg}*g^{bd}+4*g^{ad}*g^{gh}*g^{ef}*g^{bc}+4*g^{bh}*g^{af}*g^{cd}*g^{eg}+4*g^{eh}*g^{ab}*g^{cg}*g^{df}-4*g^{dh}*g^{ag}*g^{ce}*g^{bf}+4*g^{ah}*g^{dg}*g^{ce}*g^{bf}-4*g^{af}*g^{dh}*g^{cg}*g^{be}-4*g^{af}*g^{gh}*g^{bd}*g^{ce}+4*g^{ab}*g^{gh}*g^{cf}*g^{de}-4*g^{ah}*g^{cg}*g^{de}*g^{bf}-4*g^{ch}*g^{af}*g^{bg}*g^{de}+4*g^{fg}*g^{bc}*g^{ah}*g^{de}-4*g^{ag}*g^{eh}*g^{bd}*g^{cf}-4*g^{ah}*g^{df}*g^{bg}*g^{ce}-4*g^{bh}*g^{ad}*g^{eg}*g^{cf}-4*g^{ag}*g^{fh}*g^{cd}*g^{be}-4*g^{dh}*g^{ab}*g^{fg}*g^{ce}+4*g^{cd}*g^{ef}*g^{ah}*g^{bg}-4*g^{af}*g^{eh}*g^{cd}*g^{bg}-4*g^{fh}*g^{ad}*g^{bc}*g^{eg}-4*g^{dh}*g^{ag}*g^{ef}*g^{bc}+4*g^{ef}*g^{bc}*g^{ah}*g^{dg}-4*g^{ch}*g^{ab}*g^{eg}*g^{df}-4*g^{af}*g^{eh}*g^{bc}*g^{dg}-4*g^{ae}*g^{fh}*g^{bd}*g^{cg}-4*g^{ae}*g^{ch}*g^{dg}*g^{bf}+4*g^{eh}*g^{ad}*g^{cf}*g^{bg}+4*g^{ae}*g^{ch}*g^{df}*g^{bg}-4*g^{dh}*g^{ab}*g^{ef}*g^{cg}-4*g^{cd}*g^{eg}*g^{ah}*g^{bf}+4*g^{af}*g^{dh}*g^{bg}*g^{ce}-4*g^{ch}*g^{ad}*g^{fg}*g^{be}+4*g^{bh}*g^{ad}*g^{fg}*g^{ce}-4*g^{bh}*g^{ac}*g^{fg}*g^{de}-4*g^{bc}*g^{eg}*g^{ah}*g^{df}-4*g^{bh}*g^{ag}*g^{cf}*g^{de}-4*g^{bh}*g^{ac}*g^{ef}*g^{dg}+4*g^{af}*g^{gh}*g^{cd}*g^{be}-4*g^{ag}*g^{fh}*g^{bc}*g^{de}+4*g^{bh}*g^{ad}*g^{ef}*g^{cg}+4*g^{ae}*g^{dh}*g^{cg}*g^{bf}-4*g^{ae}*g^{dh}*g^{fg}*g^{bc}+4*g^{ab}*g^{fh}*g^{dg}*g^{ce}+4*g^{bd}*g^{eg}*g^{cf}*g^{ah}-4*g^{ad}*g^{gh}*g^{cf}*g^{be}-4*g^{ab}*g^{fh}*g^{cg}*g^{de}+4*g^{dh}*g^{ac}*g^{fg}*g^{be}+4*g^{ch}*g^{ag}*g^{bd}*g^{ef}-4*g^{ab}*g^{gh}*g^{df}*g^{ce}+4*g^{af}*g^{dh}*g^{bc}*g^{eg}+4*g^{bh}*g^{ac}*g^{eg}*g^{df}+4*g^{ab}*g^{gh}*g^{cd}*g^{ef}+4*g^{dh}*g^{ag}*g^{cf}*g^{be}-4*g^{cf}*g^{ah}*g^{dg}*g^{be}+4*g^{eh}*g^{ac}*g^{dg}*g^{bf}-4*g^{ac}*g^{gh}*g^{bd}*g^{ef}-4*g^{ch}*g^{ad}*g^{ef}*g^{bg}-4*g^{eh}*g^{ad}*g^{cg}*g^{bf}+4*g^{eh}*g^{ad}*g^{fg}*g^{bc}-4*g^{eh}*g^{ac}*g^{df}*g^{bg}-4*g^{fg}*g^{bd}*g^{ah}*g^{ce}-4*g^{ab}*g^{fh}*g^{cd}*g^{eg}+4*g^{af}*g^{gh}*g^{bc}*g^{de}+4*g^{ae}*g^{bh}*g^{cf}*g^{dg}+4*g^{fh}*g^{ac}*g^{bd}*g^{eg}+4*g^{ae}*g^{fh}*g^{cd}*g^{bg}+4*g^{ae}*g^{fh}*g^{bc}*g^{dg}-4*g^{ae}*g^{gh}*g^{cd}*g^{bf}-4*g^{eh}*g^{ab}*g^{cf}*g^{dg}-4*g^{bd}*g^{ef}*g^{ah}*g^{cg}+4*g^{ch}*g^{ad}*g^{eg}*g^{bf}+4*g^{af}*g^{eh}*g^{bd}*g^{cg}+4*g^{ae}*g^{ch}*g^{fg}*g^{bd}-4*g^{ch}*g^{ag}*g^{be}*g^{df}+4*g^{bh}*g^{ag}*g^{df}*g^{ce}-4*g^{ae}*g^{gh}*g^{bc}*g^{df}+4*g^{cf}*g^{ah}*g^{bg}*g^{de}-4*g^{fh}*g^{ac}*g^{dg}*g^{be}+4*g^{dh}*g^{ac}*g^{ef}*g^{bg}-4*g^{bh}*g^{ag}*g^{cd}*g^{ef}+4*g^{fh}*g^{ad}*g^{cg}*g^{be}-4*g^{ch}*g^{af}*g^{bd}*g^{eg}+4*g^{ac}*g^{gh}*g^{be}*g^{df}+4*g^{ae}*g^{gh}*g^{bd}*g^{cf}-4*g^{dh}*g^{ac}*g^{eg}*g^{bf}+4*g^{ag}*g^{fh}*g^{bd}*g^{ce}+4*g^{fg}*g^{cd}*g^{ah}*g^{be}+4*g^{ch}*g^{af}*g^{dg}*g^{be}-4*g^{bh}*g^{af}*g^{dg}*g^{ce}+4*g^{ch}*g^{ag}*g^{de}*g^{bf}+4*g^{ah}*g^{cg}*g^{be}*g^{df}+4*g^{bh}*g^{af}*g^{cg}*g^{de}+4*g^{ag}*g^{eh}*g^{cd}*g^{bf}-4*g^{fh}*g^{ad}*g^{bg}*g^{ce}+4*g^{fh}*g^{ac}*g^{bg}*g^{de}+4*g^{ch}*g^{ab}*g^{fg}*g^{de}+4*g^{dh}*g^{ab}*g^{eg}*g^{cf}-4*g^{ae}*g^{bh}*g^{fg}*g^{cd}+4*g^{ad}*g^{gh}*g^{ce}*g^{bf}-4*g^{ac}*g^{gh}*g^{de}*g^{bf}+4*g^{ag}*g^{eh}*g^{bc}*g^{df}+4*g^{eh}*g^{ab}*g^{fg}*g^{cd}-4*g^{ae}*g^{dh}*g^{cf}*g^{bg}+4*g^{ch}*g^{ab}*g^{ef}*g^{dg}-4*g^{ae}*g^{bh}*g^{cg}*g^{df}");
    }

    @Test
    public void test5() {
        GeneralIndicesInsertion indicesInsertion = new GeneralIndicesInsertion();
        CC.current().getParseManager().defaultParserPreprocessors.add(indicesInsertion);
        indicesInsertion.addInsertionRule(parseSimple("G^a'_b'a"), IndexType.LatinLower1);

        Tensor t = parse("Tr[G^a*G^b*G^c*G^d*G^e*G^f*G^g*G^h*G^i*G^j]");
        t = Expand.expand(DiracTrace.trace(t));
        Assert.assertEquals(t.size(), 945);
    }

}
