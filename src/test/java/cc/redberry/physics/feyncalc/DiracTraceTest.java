package cc.redberry.physics.feyncalc;

import cc.redberry.core.TAssert;
import cc.redberry.core.context.CC;
import cc.redberry.core.indices.IndexType;
import cc.redberry.core.parser.preprocessor.GeneralIndicesInsertion;
import cc.redberry.core.tensor.Expression;
import cc.redberry.core.tensor.SimpleTensor;
import cc.redberry.core.tensor.Tensor;
import cc.redberry.core.tensor.Tensors;
import cc.redberry.core.tensor.iterator.TensorLastIterator;
import cc.redberry.core.transformations.ContractIndices;
import cc.redberry.core.transformations.expand.Expand;
import junit.framework.Assert;
import org.junit.Test;

import static cc.redberry.core.tensor.Tensors.*;
import static cc.redberry.physics.feyncalc.DiracTrace.trace;

/**
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public class DiracTraceTest {
    @Test
    public void test1() {
        for (int i = 0; i < 100; ++i) {
            CC.resetTensorNames();
            GeneralIndicesInsertion indicesInsertion = new GeneralIndicesInsertion();
            CC.current().getParseManager().defaultParserPreprocessors.add(indicesInsertion);
            indicesInsertion.addInsertionRule(parseSimple("G^a'_b'a"), IndexType.LatinLower1);

            Tensor t = parse("Tr[G_a*G_b*G_c*G_d]");
            t = Expand.expand(trace(t));
            Tensor expected = parse("-4*g_{ac}*g_{bd}+4*g_{ad}*g_{bc}+4*g_{ab}*g_{cd}");
            TAssert.assertEquals(t, expected);
        }
    }

    @Test
    public void test2() {
        for (int i = 0; i < 100; ++i) {
            CC.resetTensorNames();
            GeneralIndicesInsertion indicesInsertion = new GeneralIndicesInsertion();
            CC.current().getParseManager().defaultParserPreprocessors.add(indicesInsertion);
            indicesInsertion.addInsertionRule(parseSimple("G^a'_b'a"), IndexType.LatinLower1);

            Tensor t = parse("Tr[G^a*G^b*G^c*G^d*G^e*G^f]");
            t = Expand.expand(trace(t));
            Tensor expected = parse("4*g^{af}*g^{be}*g^{cd}-4*g^{ae}*g^{bf}*g^{cd}+4*g^{ab}*g^{cd}*g^{ef}-4*g^{af}*g^{bd}*g^{ce}+4*g^{ad}*g^{bf}*g^{ce}+4*g^{ae}*g^{bd}*g^{cf}-4*g^{ad}*g^{be}*g^{cf}+4*g^{af}*g^{bc}*g^{de}-4*g^{ac}*g^{bf}*g^{de}+4*g^{ab}*g^{cf}*g^{de}-4*g^{ae}*g^{bc}*g^{df}+4*g^{ac}*g^{be}*g^{df}-4*g^{ab}*g^{ce}*g^{df}+4*g^{ad}*g^{bc}*g^{ef}-4*g^{ac}*g^{bd}*g^{ef}");
            TAssert.assertEquals(t, expected);
        }
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
        t = Expand.expand(trace(t));
        TAssert.assertEquals(t, "-4*g^{eh}*g^{ac}*g^{fg}*g^{bd}+4*g^{ad}*g^{gh}*g^{ef}*g^{bc}+4*g^{bh}*g^{af}*g^{cd}*g^{eg}+4*g^{eh}*g^{ab}*g^{cg}*g^{df}-4*g^{dh}*g^{ag}*g^{ce}*g^{bf}+4*g^{ah}*g^{dg}*g^{ce}*g^{bf}-4*g^{af}*g^{dh}*g^{cg}*g^{be}-4*g^{af}*g^{gh}*g^{bd}*g^{ce}+4*g^{ab}*g^{gh}*g^{cf}*g^{de}-4*g^{ah}*g^{cg}*g^{de}*g^{bf}-4*g^{ch}*g^{af}*g^{bg}*g^{de}+4*g^{fg}*g^{bc}*g^{ah}*g^{de}-4*g^{ag}*g^{eh}*g^{bd}*g^{cf}-4*g^{ah}*g^{df}*g^{bg}*g^{ce}-4*g^{bh}*g^{ad}*g^{eg}*g^{cf}-4*g^{ag}*g^{fh}*g^{cd}*g^{be}-4*g^{dh}*g^{ab}*g^{fg}*g^{ce}+4*g^{cd}*g^{ef}*g^{ah}*g^{bg}-4*g^{af}*g^{eh}*g^{cd}*g^{bg}-4*g^{fh}*g^{ad}*g^{bc}*g^{eg}-4*g^{dh}*g^{ag}*g^{ef}*g^{bc}+4*g^{ef}*g^{bc}*g^{ah}*g^{dg}-4*g^{ch}*g^{ab}*g^{eg}*g^{df}-4*g^{af}*g^{eh}*g^{bc}*g^{dg}-4*g^{ae}*g^{fh}*g^{bd}*g^{cg}-4*g^{ae}*g^{ch}*g^{dg}*g^{bf}+4*g^{eh}*g^{ad}*g^{cf}*g^{bg}+4*g^{ae}*g^{ch}*g^{df}*g^{bg}-4*g^{dh}*g^{ab}*g^{ef}*g^{cg}-4*g^{cd}*g^{eg}*g^{ah}*g^{bf}+4*g^{af}*g^{dh}*g^{bg}*g^{ce}-4*g^{ch}*g^{ad}*g^{fg}*g^{be}+4*g^{bh}*g^{ad}*g^{fg}*g^{ce}-4*g^{bh}*g^{ac}*g^{fg}*g^{de}-4*g^{bc}*g^{eg}*g^{ah}*g^{df}-4*g^{bh}*g^{ag}*g^{cf}*g^{de}-4*g^{bh}*g^{ac}*g^{ef}*g^{dg}+4*g^{af}*g^{gh}*g^{cd}*g^{be}-4*g^{ag}*g^{fh}*g^{bc}*g^{de}+4*g^{bh}*g^{ad}*g^{ef}*g^{cg}+4*g^{ae}*g^{dh}*g^{cg}*g^{bf}-4*g^{ae}*g^{dh}*g^{fg}*g^{bc}+4*g^{ab}*g^{fh}*g^{dg}*g^{ce}+4*g^{bd}*g^{eg}*g^{cf}*g^{ah}-4*g^{ad}*g^{gh}*g^{cf}*g^{be}-4*g^{ab}*g^{fh}*g^{cg}*g^{de}+4*g^{dh}*g^{ac}*g^{fg}*g^{be}+4*g^{ch}*g^{ag}*g^{bd}*g^{ef}-4*g^{ab}*g^{gh}*g^{df}*g^{ce}+4*g^{af}*g^{dh}*g^{bc}*g^{eg}+4*g^{bh}*g^{ac}*g^{eg}*g^{df}+4*g^{ab}*g^{gh}*g^{cd}*g^{ef}+4*g^{dh}*g^{ag}*g^{cf}*g^{be}-4*g^{cf}*g^{ah}*g^{dg}*g^{be}+4*g^{eh}*g^{ac}*g^{dg}*g^{bf}-4*g^{ac}*g^{gh}*g^{bd}*g^{ef}-4*g^{ch}*g^{ad}*g^{ef}*g^{bg}-4*g^{eh}*g^{ad}*g^{cg}*g^{bf}+4*g^{eh}*g^{ad}*g^{fg}*g^{bc}-4*g^{eh}*g^{ac}*g^{df}*g^{bg}-4*g^{fg}*g^{bd}*g^{ah}*g^{ce}-4*g^{ab}*g^{fh}*g^{cd}*g^{eg}+4*g^{af}*g^{gh}*g^{bc}*g^{de}+4*g^{ae}*g^{bh}*g^{cf}*g^{dg}+4*g^{fh}*g^{ac}*g^{bd}*g^{eg}+4*g^{ae}*g^{fh}*g^{cd}*g^{bg}+4*g^{ae}*g^{fh}*g^{bc}*g^{dg}-4*g^{ae}*g^{gh}*g^{cd}*g^{bf}-4*g^{eh}*g^{ab}*g^{cf}*g^{dg}-4*g^{bd}*g^{ef}*g^{ah}*g^{cg}+4*g^{ch}*g^{ad}*g^{eg}*g^{bf}+4*g^{af}*g^{eh}*g^{bd}*g^{cg}+4*g^{ae}*g^{ch}*g^{fg}*g^{bd}-4*g^{ch}*g^{ag}*g^{be}*g^{df}+4*g^{bh}*g^{ag}*g^{df}*g^{ce}-4*g^{ae}*g^{gh}*g^{bc}*g^{df}+4*g^{cf}*g^{ah}*g^{bg}*g^{de}-4*g^{fh}*g^{ac}*g^{dg}*g^{be}+4*g^{dh}*g^{ac}*g^{ef}*g^{bg}-4*g^{bh}*g^{ag}*g^{cd}*g^{ef}+4*g^{fh}*g^{ad}*g^{cg}*g^{be}-4*g^{ch}*g^{af}*g^{bd}*g^{eg}+4*g^{ac}*g^{gh}*g^{be}*g^{df}+4*g^{ae}*g^{gh}*g^{bd}*g^{cf}-4*g^{dh}*g^{ac}*g^{eg}*g^{bf}+4*g^{ag}*g^{fh}*g^{bd}*g^{ce}+4*g^{fg}*g^{cd}*g^{ah}*g^{be}+4*g^{ch}*g^{af}*g^{dg}*g^{be}-4*g^{bh}*g^{af}*g^{dg}*g^{ce}+4*g^{ch}*g^{ag}*g^{de}*g^{bf}+4*g^{ah}*g^{cg}*g^{be}*g^{df}+4*g^{bh}*g^{af}*g^{cg}*g^{de}+4*g^{ag}*g^{eh}*g^{cd}*g^{bf}-4*g^{fh}*g^{ad}*g^{bg}*g^{ce}+4*g^{fh}*g^{ac}*g^{bg}*g^{de}+4*g^{ch}*g^{ab}*g^{fg}*g^{de}+4*g^{dh}*g^{ab}*g^{eg}*g^{cf}-4*g^{ae}*g^{bh}*g^{fg}*g^{cd}+4*g^{ad}*g^{gh}*g^{ce}*g^{bf}-4*g^{ac}*g^{gh}*g^{de}*g^{bf}+4*g^{ag}*g^{eh}*g^{bc}*g^{df}+4*g^{eh}*g^{ab}*g^{fg}*g^{cd}-4*g^{ae}*g^{dh}*g^{cf}*g^{bg}+4*g^{ch}*g^{ab}*g^{ef}*g^{dg}-4*g^{ae}*g^{bh}*g^{cg}*g^{df}");
    }

    @Test
    public void test5() {
        GeneralIndicesInsertion indicesInsertion = new GeneralIndicesInsertion();
        CC.current().getParseManager().defaultParserPreprocessors.add(indicesInsertion);
        indicesInsertion.addInsertionRule(parseSimple("G^a'_b'a"), IndexType.LatinLower1);

        Tensor t = parse("Tr[G^a*G^b*G^c*G^d*G^e*G^f*G^g*G^h*G^i*G^j]");
        t = Expand.expand(trace(t));
        Assert.assertEquals(t.size(), 945);
    }

    @Test
    public void test6() {
        GeneralIndicesInsertion indicesInsertion = new GeneralIndicesInsertion();
        CC.current().getParseManager().defaultParserPreprocessors.add(indicesInsertion);
        indicesInsertion.addInsertionRule(parseSimple("G^a'_b'a"), IndexType.LatinLower1);

        Tensor t = parse("Tr[G_a*G_b*G^b*G_d] + Tr[G_a*G_d]");
        t = Expand.expand(trace(t));
        t = parseExpression("d^a_a = 4").transform(t);
        Tensor expected = parse("20*g_ad");
        TAssert.assertEquals(t, expected);
    }

    @Test
    public void test7() {
        Tensor t = parse("(4*M**(-4)*N*m+M**(-3)*N+(1/8)*((1/16)*M**4+m**4+(1/2)*M**2*m**2)**(-1)*M*N)*G^{ge'}_{a'}*G_{a}^{b'}_{d'}*G_{c}^{d'}_{e'}*G^{fh'}_{i'}*G_{m}^{i'}_{b'}*G_{n}^{a'}_{h'}*k2_{f}*k2^{a}*k2^{c}*k2_{g}*g_{AB}*k1^{m}");
        t = trace(t);
        assertContainsGamma(t);
    }

    @Test
    public void test8() {
        Tensor t = parse("g_{AB}*d^{A'}_{A'}*k2_{b}*k2_{f}*P^{b}*P^{a}*P_{g}*P_{d}*G^{fh'}_{i'}*G^{db'}_{c'}*G_{a}^{c'}_{d'}*G^{ge'}_{a'}*G_{m}^{i'}_{b'}*G_{n}^{a'}_{h'}*d^{d'}_{e'}");
        t = trace(t);
        assertContainsGamma(t);
    }

    @Test
    public void test9() {
        GeneralIndicesInsertion indicesInsertion = new GeneralIndicesInsertion();
        CC.current().getParseManager().defaultParserPreprocessors.add(indicesInsertion);
        indicesInsertion.addInsertionRule(parseSimple("G^a'_b'a"), IndexType.LatinLower1);
        Tensor t;
        t = parse("Tr[p2_a*G^a*e2_b*G^b*p1_c*G^c*e1_d*G^d*p1_e*G^e*e2_f*G^f*p1_g*G^g*e1_h*G^h]");
        TAssert.assertEquals(tr(t), "-32*p1_{g}*p1^{g}*p1^{a}*p2_{a}");
        t = parse("Tr[p2_a*G^a*e2_b*G^b*p1_c*G^c*e1_d*G^d*p1_e*G^e*e2_f*G^f*k1_g*G^g*e1_h*G^h]");
        TAssert.assertEquals(tr(t), "-32*p1_{g}*k1^{g}*p1^{a}*p2_{a}");
        t = parse("Tr[p2_a*G^a*e2_b*G^b*k1_c*G^c*e1_d*G^d*p1_e*G^e*e2_f*G^f*p1_g*G^g*e1_h*G^h]");
        TAssert.assertEquals(tr(t), "-32*p1_{g}*k1^{g}*p1^{a}*p2_{a}");
        t = parse("Tr[p2_a*G^a*e2_b*G^b*k1_c*G^c*e1_d*G^d*p1_e*G^e*e2_f*G^f*k1_g*G^g*e1_h*G^h]");
        TAssert.assertEquals(parseExpression("k1_a*k1^a = 0").transform(tr(t)), "0");
        t = parse("Tr[p2_a*G^a*e2_b*G^b*k1_c*G^c*e1_d*G^d*p1_e*G^e*e1_f*G^f*k2_g*G^g*e2_h*G^h]");
        TAssert.assertEquals(tr(t), tr("4*Tr[p2_a*G^a*k2_b*G^b*p1_c*G^c*k1_d*G^d]"));
        t = parse("Tr[p2_a*G^a*e2_b*G^b*p1_c*G^c*e1_d*G^d*p1_e*G^e*e1_f*G^f*k2_g*G^g*e2_h*G^h]");
        TAssert.assertEquals(parseExpression("p1_a*p1^a = 0").transform(tr(t)), "0");
        t = parse("Tr[p2_a*G^a*e2_b*G^b*k1_c*G^c*e1_d*G^d*p1_e*G^e*e1_f*G^f*p1_g*G^g*e2_h*G^h]");
        TAssert.assertEquals(parseExpression("p1_a*p1^a = 0").transform(tr(t)), "0");
        t = parse("Tr[p2_a*G^a*e1_b*G^b*p1_c*G^c*e2_d*G^d*p1_e*G^e*e2_f*G^f*p1_g*G^g*e1_h*G^h]");
        TAssert.assertEquals(parseExpression("p1_a*p1^a = 0").transform(tr(t)), "0");
        t = parse("Tr[p2_a*G^a*e1_b*G^b*p1_c*G^c*e2_d*G^d*p1_e*G^e*e2_f*G^f*k1_g*G^g*e1_h*G^h]");
        TAssert.assertEquals(parseExpression("p1_a*p1^a = 0").transform(tr(t)), "0");
        t = parse("Tr[p2_a*G^a*e1_b*G^b*k2_c*G^c*e2_d*G^d*p1_e*G^e*e2_f*G^f*p1_g*G^g*e1_h*G^h]");
        TAssert.assertEquals(parseExpression("p1_a*p1^a = 0").transform(tr(t)), "0");
        t = parse("Tr[p2_a*G^a*e1_b*G^b*k2_c*G^c*e2_d*G^d*p1_e*G^e*e2_f*G^f*k1_g*G^g*e1_h*G^h]");
        TAssert.assertEquals(tr(t), tr("4*Tr[p2_a*G^a*k2_b*G^b*p1_c*G^c*k1_d*G^d]"));
        t = parse("Tr[p2_a*G^a*e1_b*G^b*p1_c*G^c*e2_d*G^d*p1_e*G^e*e1_f*G^f*p1_g*G^g*e2_h*G^h]");
        TAssert.assertEquals(parseExpression("p1_a*p1^a = 0").transform(tr(t)), "0");
        t = parse("Tr[p2_a*G^a*e1_b*G^b*p1_c*G^c*e2_d*G^d*p1_e*G^e*e1_f*G^f*k2_g*G^g*e2_h*G^h]");
        TAssert.assertEquals(parseExpression("p1_a*p1^a = 0").transform(tr(t)), "-32*p1_{g}*k2^{g}*p1^{a}*p2_{a}");
        t = parse("Tr[p2_a*G^a*e1_b*G^b*k2_c*G^c*e2_d*G^d*p1_e*G^e*e1_f*G^f*k2_g*G^g*e2_h*G^h]");
        TAssert.assertEquals(parseExpression("k2_a*k2^a = 0").transform(tr(t)), "0");
    }

    private static Tensor tr(String t) {
        return tr(parse(t));
    }

    private static Tensor tr(Tensor t) {
        t = trace(t);
        t = parseExpression("e1_m*e1_n = g_mn").transform(t);
        t = parseExpression("e2_m*e2_n = g_mn").transform(t);
        t = ContractIndices.contract(t);
        t = parseExpression("d^m_m = 4").transform(t);
        return t;
    }

    @Test
    public void test10() {
        GeneralIndicesInsertion indicesInsertion = new GeneralIndicesInsertion();
        CC.current().getParseManager().defaultParserPreprocessors.add(indicesInsertion);
        indicesInsertion.addInsertionRule(parseSimple("G^a'_b'a"), IndexType.LatinLower1);
        indicesInsertion.addInsertionRule(parseSimple("G5^a'_b'"), IndexType.LatinLower1);
        Tensor t;
        t = parse("Tr[G_a*G5]");
        TAssert.assertEquals(trace(t), "0");
        t = parse("Tr[G_a*G_b*G5]");
        TAssert.assertEquals(trace(t), "0");
        t = parse("Tr[G_a*G_b*G_c*G5]");
        TAssert.assertEquals(trace(t), "0");
        t = parse("Tr[G_a*G_b*G_c*G_d*G_e*G5]");
        TAssert.assertEquals(trace(t), "0");
        t = parse("Tr[G_a*G_b*G_c*G_d*G_e*G_h*G_f*G5]");
        TAssert.assertEquals(trace(t), "0");
    }

    @Test
    public void test11() {
        GeneralIndicesInsertion indicesInsertion = new GeneralIndicesInsertion();
        CC.current().getParseManager().defaultParserPreprocessors.add(indicesInsertion);
        indicesInsertion.addInsertionRule(parseSimple("G^a'_b'a"), IndexType.LatinLower1);
        indicesInsertion.addInsertionRule(parseSimple("G5^a'_b'"), IndexType.LatinLower1);
        setAntiSymmetric(parseSimple("e_abcd"));
        Tensor t;
        t = parse("Tr[G_a*G_b*G_c*G_d*G5]");
        TAssert.assertEquals(trace(t), "(-4*I)*e_{abcd}");
        t = parse("Tr[G5*G_a*G_b*G5*G_c*G_d*G5]");
        TAssert.assertEquals(trace(t), "(-4*I)*e_{abcd}");
        t = parse("Tr[G5*G_a*G_b*G5*G_c*G_d*G5*G5*G5]");
        TAssert.assertEquals(trace(t), "(-4*I)*e_{abcd}");
        t = parse("Tr[G5*G_a*G5*G5*G_b*G5*G_c*G5*G5*G_d*G5*G5*G5]");
        TAssert.assertEquals(trace(t), "(-4*I)*e_{abcd}");
        t = parse("Tr[G5*G5*G_a*G5*G5*G_b*G5*G_c*G5*G5*G_d*G5*G5*G5]");
        TAssert.assertEquals(trace(t), trace(parse("Tr[G_a*G_b*G_c*G_d]")));
        t = parse("a*Tr[G5*G5*G_a*G5*G5*G_b*G5*G_c*G5*G5*b*G_d*G5*G5*G5]");
        TAssert.assertEquals(trace(t), trace(parse("a*b*Tr[G_a*G_b*G_c*G_d]")));
        t = parse("Tr[G5*G5*G5*G5*G5*G5]");
        TAssert.assertEquals(trace(t), "4");
        t = parse("Tr[G5*G5*G5*G5*G5]");
        TAssert.assertEquals(trace(t), "0");
    }


    @Test
    public void test12() {
        GeneralIndicesInsertion indicesInsertion = new GeneralIndicesInsertion();
        CC.current().getParseManager().defaultParserPreprocessors.add(indicesInsertion);
        indicesInsertion.addInsertionRule(parseSimple("G^a'_b'a"), IndexType.LatinLower1);
        indicesInsertion.addInsertionRule(parseSimple("G5^a'_b'"), IndexType.LatinLower1);
        setAntiSymmetric(parseSimple("e_abcd"));
        Tensor t;

        t = parse("Tr[G_a*G_b*G_c*G_d*G_e*G^d*G5]");
        TAssert.assertEquals(trace(t), "8*I*e_abce");
        t = parse("Tr[G_a*G_b*G_c*G_d*G_e*G^a*G5]");
        TAssert.assertEquals(trace(t), "16*I*e_bcde");
        t = parse("Tr[G_a*G_b*G_c*G_d*G_e*G^b*G5]");
        TAssert.assertEquals(trace(t), "-8*I*e_acde");
        t = parse("Tr[G_a*G_b*G_c*G_d*G_e*G^c*G5]");
        TAssert.assertEquals(trace(t), "0");
        t = parse("Tr[G_a*G_b*G_c*G_d*G_e*G^e*G5]");
        TAssert.assertEquals(trace(t), "-16*I*e_abcd");
        t = parse("-Tr[G5*G_a*G5*G_b*G_c*G_d*G_e*G^d*G5]");
        TAssert.assertEquals(trace(t), "8*I*e_abce");
        t = parse("Tr[G_a*G5*G_b*G_c*G5*G5*G5*G_d*G_e*G^a*G5]");
        TAssert.assertEquals(trace(t), "16*I*e_bcde");
        t = parse("Tr[G_a*G5*G_b*G_c*G_d*G5*G_e*G^b*G5]");
        TAssert.assertEquals(trace(t), "8*I*e_acde");
        t = parse("Tr[-G5*G_a*G_b*G_c*G5*G_d*G_e*G^c*G5]");
        TAssert.assertEquals(trace(t), "0");
        t = parse("-Tr[G5*G_a*G5*G_b*G5*G_c*G_d*G_e*G^e*G5*G5]");
        TAssert.assertEquals(trace(t), "-16*I*e_abcd");
    }

    @Test
    public void test13() {
        GeneralIndicesInsertion indicesInsertion = new GeneralIndicesInsertion();
        CC.current().getParseManager().defaultParserPreprocessors.add(indicesInsertion);
        indicesInsertion.addInsertionRule(parseSimple("G^a'_b'a"), IndexType.LatinLower1);
        indicesInsertion.addInsertionRule(parseSimple("G5^a'_b'"), IndexType.LatinLower1);
        setAntiSymmetric(parseSimple("e_abcd"));
        Tensor t, expected;

        t = parse("Tr[G_a*G_b*G_c*G_d*G_e*G_f*G_g*G^a*G5]");
        expected = parse("16*I*g_{bf}*e_{cdeg}-16*I*g_{bg}*e_{cdef}+16*I*g_{cd}*e_{befg}-16*I*g_{ce}*e_{bdfg}+16*I*g_{de}*e_{bcfg}+16*I*g_{fg}*e_{bcde}");
        TAssert.assertEquals(trace(t), expected);
        t = parse("Tr[G_a*G_b*G_c*G_d*G_e*G^d*G_g*G^b*G5]");
        TAssert.assertEquals(trace(t), "16*I*e_aceg");
    }

    //Expression shouten = parseExpression("g_fa*e_bcde = -(g_fb*e_cdea + g_fc*e_deab+g_fd*e_eabc+g_fe*e_abcd)");
    private static final SimpleTensor defaultGamma = parseSimple("G^a'_b'a");

    private static void assertContainsGamma(Tensor t) {
        TensorLastIterator iterator = new TensorLastIterator(t);
        Tensor c;
        while ((c = iterator.next()) != null) {
            if (c instanceof SimpleTensor && ((SimpleTensor) c).getName() == defaultGamma.getName())
                throw new AssertionError();
        }
    }

}
