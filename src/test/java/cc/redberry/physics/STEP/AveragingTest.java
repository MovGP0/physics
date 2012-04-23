package cc.redberry.physics.STEP;

import cc.redberry.core.context.CC;
import cc.redberry.core.indices.IndexType;
import cc.redberry.core.tensor.Expression;
import cc.redberry.core.tensor.Sum;
import cc.redberry.core.tensor.Tensor;
import cc.redberry.transformation.CalculateNumbers;
import cc.redberry.transformation.ExpandBrackets;
import cc.redberry.transformation.Transformations;
import cc.redberry.transformation.Transformer;
import cc.redberry.transformation.collect.CollectFactory;
import cc.redberry.transformation.contractions.IndicesContractionsTransformation;
import org.junit.Test;

/**
 * Created by IntelliJ IDEA. User: Konstantin_2 Date: 15.04.12 Time: 22:59 To
 * change this template use File | Settings | File Templates.
 */
public class AveragingTest {
//    @Test
//    public void test1() {
//        Tensor t = Averaging.average(new int[]{0, 1, 2, 3, 4, 5, 6, 7});
//        t = Transformations.expandBrackets(t);
//        System.out.println(t);
//        System.out.println("Length: " + ((Sum) t).size());
//        t = CollectFactory.createCollectAllEqualTerms().transform(t);
//        System.out.println(t);
//        System.out.println("Length 2: " + ((Sum) t).size());
//    }
    @Test
    public void test2() {
        Tensor t = CC.parse("n_{\\nu}*n_{\\alpha}*n^{\\beta}*n^{\\gamma}");
        t = Averaging.INSTANCE.transform(t);
        System.out.println(t);
    }

    @Test
    public void test3() {
        Tensor t = CC.parse("a*n_\\mu");
        t = Averaging.INSTANCE.transform(t);
        System.out.println(t);
    }

    @Test
    public void test4() {
        Tensor t = CC.parse("n^\\mu*n_\\mu*n_\\alpha*n^\\alpha*n_\\nu*n^\\nu*n_\\lambda*n^\\lambda*n_\\rho*n^\\rho");
        Expression d = new Expression("d_\\mu^\\mu=4");
        t = Averaging.INSTANCE.transform(t);
        t = IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC.transform(t);
        t = d.transform(t);
        t = CalculateNumbers.INSTANCE.transform(t);
        System.out.println(t);
    }

    @Test
    public void test5() {
        Expression ff = new Expression("FF=(-1/6)*F^{\\nu \\beta \\epsilon }_{\\zeta }*F_{\\nu \\beta }^{\\zeta }_{\\epsilon }+n^{\\mu }*F^{\\alpha }_{\\nu }^{\\epsilon }_{\\lambda }*n^{\\nu }*F_{\\alpha \\mu }^{\\lambda }_{\\epsilon }+(-8/3)*n^{\\mu }*F_{\\beta \\nu }^{\\epsilon }_{\\lambda }*n^{\\alpha }*n^{\\beta }*n^{\\nu }*F_{\\alpha \\mu }^{\\lambda }_{\\epsilon }");
        CC.addSymmetry("F_{\\mu\\nu\\alpha\\beta}", IndexType.GreekLower, true, new int[]{1, 0, 2, 3});
        ff.eval(
                Averaging.INSTANCE,
                new Transformer(ExpandBrackets.EXPAND_ALL),
                IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                new Expression("F_{\\mu}^\\mu_\\alpha\\beta=0"),
                CalculateNumbers.INSTANCE,
                CollectFactory.createCollectAllEqualTerms(),
                CalculateNumbers.INSTANCE);
        System.out.println(ff);
    }
}
