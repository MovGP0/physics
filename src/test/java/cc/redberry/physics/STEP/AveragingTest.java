package cc.redberry.physics.STEP;

import cc.redberry.core.context.CC;
import cc.redberry.core.tensor.Sum;
import cc.redberry.core.tensor.Tensor;
import cc.redberry.transformation.Transformations;
import cc.redberry.transformation.collect.CollectFactory;
import org.junit.Test;

/**
 * Created by IntelliJ IDEA.
 * User: Konstantin_2
 * Date: 15.04.12
 * Time: 22:59
 * To change this template use File | Settings | File Templates.
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
        Tensor t = CC.parse("a_a*a^a");
        t = Averaging.INSTANCE.transform(t);
        System.out.println(t);
    }
}
