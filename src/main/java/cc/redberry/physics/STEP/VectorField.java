package cc.redberry.physics.STEP;

import cc.redberry.core.context.CC;
import cc.redberry.core.tensor.Expression;
import cc.redberry.core.tensor.Tensor;
import cc.redberry.core.tensor.testing.TTest;
import cc.redberry.core.transformations.RenameConflictingIndices;
import cc.redberry.core.utils.Indicator;
import cc.redberry.transformation.*;
import cc.redberry.transformation.collect.CollectFactory;
import cc.redberry.transformation.collect.CollectPowers;
import cc.redberry.transformation.concurrent.EACScalars;
import cc.redberry.transformation.contractions.IndicesContractionsTransformation;
import cc.redberry.transformation.substitutions.TensorTreeIndicatorImpl;

import static cc.redberry.core.indices.IndexType.GreekLower;
import static cc.redberry.physics.util.IndicesFactoryUtil.createIndices;

/**
 * Created by IntelliJ IDEA.
 * User: Konstantin_2
 * Date: 08.04.12
 * Time: 12:36
 * To change this template use File | Settings | File Templates.
 */
public class VectorField extends MainTensors {
    public final Expression L = new Expression("L = 2");

    public final Expression KRONECKER_DIMENSION =
            new Expression("d^{\\alpha}_{\\alpha} = 4");

    public final Expression Kn =
            new Expression("Kn_\\alpha^\\beta=d_\\alpha^\\beta-\\lambda*n_\\alpha*n^\\beta");
    public final Expression KINV =
            new Expression("KINV_\\alpha^\\beta=d_\\alpha^\\beta+\\gamma*n_\\alpha*n^\\beta");
    public static final Tensor[] MATRIX_Kn = {
            CC.parse("Kn_\\alpha^\\beta"),
            CC.parse("KINV_\\alpha^\\beta")
    };
    //public final Expression Gamma = new Expression("\\gamma=\\frac*\\lambda*1-\\lambda");
    public final Expression Gamma = new Expression("\\gamma=0");
    public final Expression Lambda = new Expression("\\lambda=0");

    public static final Expression K_2 =
            new Expression("K^{\\mu\\nu}_\\alpha^\\beta=g^{\\mu\\nu}*d_{\\alpha}^{\\beta}-\\lambda/2*(g^{\\mu\\beta}*d_\\alpha^\\nu+g^{\\nu\\beta}*d_\\alpha^\\mu)");
     public static final Expression HATW = new Expression("HATW^{\\alpha}_{\\beta}=P^{\\alpha}_{\\beta}+\\lambda/2*R^{\\alpha}_{\\beta}");

    public static final Tensor[] MATRIX_K = {
            CC.parse("K^{\\mu\\nu}")};

    public static void main(String[] args) {
        VectorField vf = new VectorField();
        vf.start();
    }

    public void start() {
        for (Expression ex : ALL)
            ex.eval(new Transformer(RenameConflictingIndices.INSTANCE));
        for (Expression ex : ALL)
            ex.eval(L.asSubstitution(), CalculateNumbers.INSTANCE);
        CC.addSymmetry("P_{\\alpha\\beta}", GreekLower, true, 1,0);
        insertIndices();
        System.out.println(ACTION);
        evalHatK();
        evalHatW();
        System.out.println(HATK_1);
        System.out.println(HATK_2);
        System.out.println();
        System.out.println("---------DALTAs-----------");
        System.out.println();
        evalHatKDelta();
        for(Expression e : DELTAs)
        System.out.println(e);
        System.out.println();
        System.out.println("---------TERMs-----------");
        System.out.println();
        evalTerms();
        System.out.println(Flat);
        System.out.println();
        System.out.println("----------ACTION----------");
        System.out.println();
        evalAction();
        System.out.println(ACTION);
        // evalRR();
        //for(Expression e : HATKs)
        //System.out.println(e);
    }

    public boolean indicesInserted = false;

    public void insertIndices() {
        if (indicesInserted)
            throw new IllegalAccessError("Indices are already inserted");
        Transformation indicesInsertion;
        indicesInsertion = new IndicesInsertion(matricesIndicator, createIndices(HATKs, "^{\\mu}_{\\nu}"));
        for (Expression hatK : HATKs)
            hatK.eval(indicesInsertion);
        indicesInsertion = new IndicesInsertion(matricesIndicator, createIndices(DELTAs, "^{\\mu}_{\\nu}"));
        for (Expression delta : DELTAs)
            delta.eval(indicesInsertion);
        indicesInsertion = new IndicesInsertion(matricesIndicator, createIndices(TERMs, "^{\\mu}_{\\nu}"));
        for (Expression term : TERMs)
            term.eval(indicesInsertion);
        indicesInserted = true;
    }

    public final void evalHatK() {
        for (Expression hatK : HATKs)
            hatK.eval(
                    K_2.asSubstitution(),
                    KINV.asSubstitution(),
                    Lambda.asSubstitution(), // lambda = 0
                    Gamma.asSubstitution(),  // gamma = 0
                    new Transformer(ExpandBrackets.EXPAND_EXCEPT_SYMBOLS),
                    IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                    KRONECKER_DIMENSION.asSubstitution(),
                    CollectFactory.createCollectEqualTerms(),
                    CalculateNumbers.INSTANCE,
                    EACScalars.getTransformer(),
                    CalculateNumbers.INSTANCE);
    }

    public final void evalHatW() {
            HATW.eval(
                    Lambda.asSubstitution(), // lambda = 0
                    Gamma.asSubstitution(),  // gamma = 0
                    new Transformer(ExpandBrackets.EXPAND_EXCEPT_SYMBOLS),
                    IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                    KRONECKER_DIMENSION.asSubstitution(),
                    CollectFactory.createCollectEqualTerms(),
                    CalculateNumbers.INSTANCE,
                    EACScalars.getTransformer(),
                    CalculateNumbers.INSTANCE);
    }

    public final void evalHatKDelta() {
        Transformation indicesInsertion;
        indicesInsertion = new IndicesInsertion(matricesIndicator, createIndices(DELTAs, "^{\\mu\\nu}_{\\alpha\\beta}"));
        for (Expression delta : DELTAs)
            delta.eval(
                    //indicesInsertion,
                    Lambda.asSubstitution(), // lambda = 0
                    Gamma.asSubstitution(),  // gamma = 0
                    L.asSubstitution(),
                    CalculateNumbers.INSTANCE,
                    HATK_1.asSubstitution(),
                    HATK_2.asSubstitution(),
                    HATK_3.asSubstitution(),
                    HATK_4.asSubstitution(),
                    new Transformer(ExpandBrackets.EXPAND_EXCEPT_SYMBOLS),
                    IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                    KRONECKER_DIMENSION.asSubstitution(),
                    CollectFactory.createCollectEqualTerms(),
                    CalculateNumbers.INSTANCE,
                    CollectFactory.createCollectAllScalars(),
                    CalculateNumbers.INSTANCE);
    }

    public final void evalDeltas() {
        Transformation indicesInsertion;
        indicesInsertion = new IndicesInsertion(matricesIndicator, createIndices(DELTAs, "^{\\mu\\nu}_{\\alpha\\beta}"));
        for (Expression delta : DELTAs)
            delta.eval(
                    //indicesInsertion,
                    L.asSubstitution(),
                    CalculateNumbers.INSTANCE,
                    new Transformer(RenameConflictingIndices.INSTANCE),
                    new Transformer(ExpandBrackets.EXPAND_EXCEPT_SYMBOLS),
                    //                    IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                    KRONECKER_DIMENSION.asSubstitution(),
                    //                    CollectFactory.createCollectEqualTerms1(),
                    CalculateNumbers.INSTANCE //                    ,
                    //                    CollectFactory.createCollectAllScalars(),
                    //                    CalculateNumbers.INSTANCE
            );
    }

    public final void evalTerms() {
        //Transformation indicesInsertion;
        //indicesInsertion = new IndicesInsertion(matricesIndicator, doubleAndDumpIndices(createIndices(TERMs, "^{\\mu\\nu}")));
        for(Expression e : TERMs)
        e.eval(
                //indicesInsertion,
                L.asSubstitution(),
                HATW.asSubstitution(),
                CalculateNumbers.INSTANCE,
                //RICCI.asSubstitution(),
                //RIMAN.asSubstitution(),
                new Transformer(RenameConflictingIndices.INSTANCE),
                new Transformer(ExpandBrackets.EXPAND_EXCEPT_SYMBOLS),
                IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                KRONECKER_DIMENSION.asSubstitution(),
                CalculateNumbers.INSTANCE,
                new Transformer(CollectPowers.INSTANCE),
                CollectFactory.createCollectEqualTerms());
    }

    public final void evalAction() {
        //Transformation indicesInsertion;
        //indicesInsertion = new IndicesInsertion(matricesIndicator, doubleAndDumpIndices(createIndices(TERMs, "^{\\mu\\nu}")));
        ACTION.eval(
                //indicesInsertion,
                L.asSubstitution(),
                Flat.asSubstitution(), WR.asSubstitution(), SR.asSubstitution(), SSR.asSubstitution(), FF.asSubstitution(), FR.asSubstitution(), RR.asSubstitution(),
                CalculateNumbers.INSTANCE,
                //RICCI.asSubstitution(),
                //RIMAN.asSubstitution(),
                new Transformer(RenameConflictingIndices.INSTANCE),
                new Transformer(ExpandBrackets.EXPAND_EXCEPT_SYMBOLS),
                IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                KRONECKER_DIMENSION.asSubstitution(),
                CalculateNumbers.INSTANCE,
                new Transformer(CollectPowers.INSTANCE),
                CollectFactory.createCollectEqualTerms());
    }

    public final Indicator<Tensor> matricesIndicator = new TensorTreeIndicatorImpl(new Indicator<Tensor>() {
        @Override
        public boolean is(Tensor tensor) {
            for (Tensor m : MATRICES)
                if (TTest.testEqualstensorStructure(tensor, m))
                    return true;
            for (Tensor m : MATRIX_K)
                if (TTest.testEqualstensorStructure(tensor, m))
                    return true;
            return false;
        }
    });

    public final Indicator<Tensor> simpleIndicator = new TensorTreeIndicatorImpl(new Indicator<Tensor>() {
        @Override
        public boolean is(Tensor tensor) {
            for(Tensor m : MATRIX_ACTION)
            if (TTest.testEqualstensorStructure(tensor, m))
                return true;
            return false;
        }
    });
}
