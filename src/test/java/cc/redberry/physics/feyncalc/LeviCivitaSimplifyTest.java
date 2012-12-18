package cc.redberry.physics.feyncalc;

import cc.redberry.core.TAssert;
import cc.redberry.core.context.ContextManager;
import cc.redberry.core.tensor.SimpleTensor;
import cc.redberry.core.tensor.Tensor;
import junit.framework.TestCase;
import org.junit.Test;

import static cc.redberry.core.number.Complex.ZERO;
import static cc.redberry.core.tensor.Tensors.*;
import static cc.redberry.physics.feyncalc.LeviCivitaSimplify.simplifyLeviCivita;

/**
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public class LeviCivitaSimplifyTest extends TestCase {

    @Test
    public void test1() {
        SimpleTensor eps = parseSimple("e_abcd");
        Tensor t;

        t = parse("e_abcd*k^a*k^b");
        TAssert.assertEquals(simplifyLeviCivita(t, eps), ZERO);

        t = parse("e_abcd*k^ac*k^be");
        System.out.println(simplifyLeviCivita(t, eps));
        TAssert.assertEquals(simplifyLeviCivita(t, eps), t);

        t = parse("e_abed*k^ac*k^b_c");
        System.out.println(simplifyLeviCivita(t, eps));
        TAssert.assertEquals(simplifyLeviCivita(t, eps), ZERO);

        t = parse("e_abed*g^ed");
        System.out.println(simplifyLeviCivita(t, eps));
        TAssert.assertEquals(simplifyLeviCivita(t, eps), ZERO);

        t = parse("e_abed*e^abpq*g^ed");
        System.out.println(simplifyLeviCivita(t, eps));
        TAssert.assertEquals(simplifyLeviCivita(t, eps), ZERO);

        t = parse("e_abed*e^abpq*(g^ek*g^dl+g^el*g^dk)");
        TAssert.assertEquals(simplifyLeviCivita(t, eps), ZERO);
    }

    @Test
    public void test2() {
        SimpleTensor eps = parseSimple("e_ab");
        Tensor t;

        t = parse("e_ed*e^pq");
        TAssert.assertEquals(simplifyLeviCivita(t, eps), "d^{p}_{d}*d^{q}_{e}-d^{q}_{d}*d^{p}_{e}");

        t = parse("e_ed*e^eq");
        TAssert.assertEquals(simplifyLeviCivita(t, eps), "-d_{d}^{q}");
    }


    @Test
    public void test3() {
        SimpleTensor eps = parseSimple("e_abc");
        Tensor t;

        t = parse("e_abc*e^abd");
        TAssert.assertEquals(simplifyLeviCivita(t, eps), "2*d^d_c");
        t = parse("e_abc*e^abc");
        TAssert.assertEquals(simplifyLeviCivita(t, eps), "6");
    }

    @Test
    public void test4() {
        SimpleTensor eps = parseSimple("e_abcf");
        addAntiSymmetry("e_abcd", 1, 0, 2, 3);
        addAntiSymmetry("e_abcd", 1, 2, 3, 0);
        Tensor t;
        t = parse("e_abcx");
        TAssert.assertEquals(simplifyLeviCivita(t, eps), t);

        t = parse("e_abcx*e^abcy");
        TAssert.assertEquals(simplifyLeviCivita(t, eps), "-6*d^y_x");
        t = parse("e_abcx*e^acby");
        TAssert.assertEquals(simplifyLeviCivita(t, eps), "6*d^y_x");
        t = parse("e_abcx*e^acby");
        TAssert.assertEquals(simplifyLeviCivita(t, eps), "6*d^y_x");
        t = parse("e_abcd*e^abcd");
        TAssert.assertEquals(simplifyLeviCivita(t, eps), "-24");
        t = parse("e_abce*e^pqrs*e_rs^ce");
        TAssert.assertEquals(simplifyLeviCivita(t, eps), "-4*e_{ab}^{pq}");
        t = parse("-4*I*e^{dh}_{b}^{f}*e_{g}^{b}_{ah}*e_{cdef}");
        TAssert.assertEquals(simplifyLeviCivita(t, eps), "16*I*e_aceg");
        t = parse("(4*I)*e^{h}_{d}^{fb}*e_{abch}*e_{e}^{d}_{gf}");
        TAssert.assertEquals(simplifyLeviCivita(t, eps), "16*I*e_aceg");
    }
}
