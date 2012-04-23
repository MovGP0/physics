package cc.redberry.physics.STEP;

import cc.redberry.core.context.CC;
import cc.redberry.core.context.ToStringMode;
import static cc.redberry.core.indices.IndexType.GreekLower;
import cc.redberry.core.number.ComplexElement;
import cc.redberry.core.number.RationalElement;
import cc.redberry.core.tensor.*;
import cc.redberry.core.tensor.testing.TTest;
import cc.redberry.core.transformations.RenameConflictingIndices;
import cc.redberry.core.utils.Indicator;
import static cc.redberry.physics.util.IndicesFactoryUtil.createIndices;
import static cc.redberry.physics.util.IndicesFactoryUtil.doubleAndDumpIndices;
import cc.redberry.physics.util.SqrSubs;
import cc.redberry.transformation.*;
import cc.redberry.transformation.collect.CollectFactory;
import cc.redberry.transformation.collect.CollectPowers;
import cc.redberry.transformation.concurrent.EACScalars;
import cc.redberry.transformation.contractions.IndicesContractionsTransformation;
import cc.redberry.transformation.substitutions.TensorTreeIndicatorImpl;

public class VectorField extends MainTensors {
    public final Expression L = new Expression("L=2");
//    public final Expression Kn =
//            new Expression("Kn_\\alpha=n_\\alpha");
//    public final Expression Kn_1 =
//            new Expression("Kn_\\alpha=Kn_\\alpha^\\beta*n_\\beta");
    
//    public static final Tensor[] MATRIX_Kn = {
//        CC.parse("Kn_\\alpha^\\beta"),
//        CC.parse("KINV_\\alpha^\\beta")
//    };
    //public final Expression Gamma = new Expression("\\gamma=\\frac*\\lambda*1-\\lambda");
    public static final Expression K_2 =
            new Expression("K^{\\mu\\nu}_\\alpha^\\beta=g^{\\mu\\nu}*d_{\\alpha}^{\\beta}");//-g^\\mu\\beta*d_\\alpha^\\nu-g^\\nu\\beta*d_\\alpha^\\mu");
    public final Expression KINV =
            new Expression("KINV_\\alpha^\\beta=d_\\alpha^\\beta");//-2*n_\\alpha*n^\\beta");
    public static final Expression W = new Expression("W^{\\alpha}_{\\beta}=P^{\\alpha}_{\\beta}");
    public static final Expression HATS_0 = new Expression("HATS=0");
    public static final Expression HATS_1 = new Expression("HATS_\\alpha=0");
    public static final Expression HATS_2 = new Expression("HATS_{\\alpha\\beta}=0");
    public static final Expression HATS_3 = new Expression("HATS_{\\alpha\\beta\\gamma}=0");
    public static final Expression HATS_4 = new Expression("HATS_{\\alpha\\beta\\gamma\\nu}=0");
    public static final Expression NABLAS_0 = new Expression("NABLAS=0");
    public static final Expression NABLAS_1 = new Expression("NABLAS_\\alpha=0");
    public static final Expression NABLAS_2 = new Expression("NABLAS_{\\alpha\\beta}=0");
    public static final Expression NABLAS_3 = new Expression("NABLAS_{\\alpha\\beta\\gamma}=0");
    public static final Expression NABLAS_4 = new Expression("NABLAS_{\\alpha\\beta\\gamma\\nu}=0");
    public static final Expression HATN_0 = new Expression("HATN=0");
    public static final Expression HATN_1 = new Expression("HATN_\\alpha=0");
    public static final Expression HATN_2 = new Expression("HATN_{\\alpha\\beta}=0");
    public static final Expression HATN_3 = new Expression("HATN_{\\alpha\\beta\\gamma}=0");
    public static final Expression HATN_4 = new Expression("HATN_{\\alpha\\beta\\gamma\\nu}=0");
    public static final Expression HATM_0 = new Expression("HATM=0");
    public static final Expression HATM_1 = new Expression("HATM_\\alpha=0");
    public static final Expression HATM_2 = new Expression("HATM_{\\alpha\\beta}=0");
    public static final Expression HATM_3 = new Expression("HATM_{\\alpha\\beta\\gamma}=0");
    public static final Expression HATM_4 = new Expression("HATM_{\\alpha\\beta\\gamma\\nu}=0");
    public static final Expression[] ZEROs = {HATS_0, HATS_1, HATS_2, HATS_3, HATS_4, HATN_0, HATN_1, HATN_2, HATN_3, HATN_4, HATM_0, HATM_1, HATM_2, HATM_3, HATM_4, NABLAS_0, NABLAS_1, NABLAS_2, NABLAS_3, NABLAS_4};
    public static final Tensor[] MATRIX_K = {
        CC.parse("K^{\\mu\\nu}"), CC.parse("W")};

    public static void main(String[] args) {
        VectorField vf = new VectorField();
        vf.start();
    }

    public void start() {
        CC.setDefaultPrintMode(ToStringMode.REDBERRY_SOUT);
        for (Expression ex : ALL)
            ex.eval(new Transformer(RenameConflictingIndices.INSTANCE));
        insertIndices();

        CC.addSymmetry("R_\\mu\\nu", GreekLower, false, new int[]{1, 0});
        CC.addSymmetry("R_\\mu\\nu\\alpha\\beta", GreekLower, true, new int[]{0, 1, 3, 2});
        CC.addSymmetry("R_\\mu\\nu\\alpha\\beta", GreekLower, false, new int[]{2, 3, 0, 1});
        CC.addSymmetry("F_\\mu\\nu\\alpha\\beta", GreekLower, true, new int[]{1, 0, 2, 3});
        CC.addSymmetry("P_\\alpha\\beta", GreekLower, false, new int[]{1, 0});
//        CC.addSymmetry("F_\\mu\\nu", GreekLower, true, new int[]{1, 0});

        for (Expression ex : ALL)
            ex.eval(L.asSubstitution(), CalculateNumbers.INSTANCE);
        System.out.println(FF);
        evalHatK();
        evalHatW();
        System.out.println(HATK_1);
        System.out.println(HATK_2);
        System.out.println(HATK_3);
        System.out.println(HATK_4);
        System.out.println();
        System.out.println("---------DALTAs-----------");
        System.out.println();
        evalHatKDelta();
        for (Expression e : DELTAs)
            System.out.println(e);
        System.out.println();
        System.out.println("---------TERMs-----------");
        System.out.println();
        evalTerms();
        System.out.println(Flat);
        System.out.println(WR);
        System.out.println(SR);
        System.out.println(SSR);
        System.out.println(FF);
        System.out.println(FR);
        System.out.println(RR);
        System.out.println();
        System.out.println("----------ACTION----------");
        System.out.println(ACTION);
        evalAction();
        System.out.println(ACTION.toString(ToStringMode.REDBERRY_SOUT));
        //latex.push(ACTION.toString());
        //latex.output();
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
        HATW.eval(indicesInsertion);
        indicesInsertion = new IndicesInsertion(matricesIndicator, createIndices(DELTAs, "^{\\mu}_{\\nu}"));
        for (Expression delta : DELTAs)
            delta.eval(indicesInsertion);
        indicesInsertion = new IndicesInsertion(matricesIndicator, doubleAndDumpIndices(createIndices(TERMs, "^{\\nu}")));
        for (Expression term : TERMs)
            term.eval(indicesInsertion);
        indicesInserted = true;
    }

    public final void evalHatK() {
        for (Expression hatK : HATKs)
            hatK.eval(
                    K_2.asSubstitution(),
                    KINV.asSubstitution(),
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
                W.asSubstitution(),
                K_2.asSubstitution(),
                KINV.asSubstitution(),
                new Transformer(ExpandBrackets.EXPAND_EXCEPT_SYMBOLS),
                IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                KRONECKER_DIMENSION.asSubstitution(),
                CollectFactory.createCollectEqualTerms(),
                CalculateNumbers.INSTANCE,
                EACScalars.getTransformer(),
                CalculateNumbers.INSTANCE);
    }

    private void evalHatKDelta() {
        for (Expression delta : DELTAs)
            delta.eval(
                    //indicesInsertion,
                    L.asSubstitution(),new Transformer(new POO()),
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
                    L.asSubstitution(),new Transformer(new POO()),
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
        for (Expression e : TERMs) {
            for (Expression d : DELTAs)
                e.eval(
                        d.asSubstitution());
            for (Expression d : HATKs)
                e.eval(
                        d.asSubstitution());

            for (Expression d : ZEROs)
                e.eval(
                        d.asSubstitution(),
                        CalculateNumbers.INSTANCE);
            e.eval(
                    L.asSubstitution(),new Transformer(new POO()),
                    //                    Kn_1.asSubstitution(),
                    //                    Kn.asSubstitution(),
                    HATW.asSubstitution(),
                    CalculateNumbers.INSTANCE,
                    //RICCI.asSubstitution(),
                    //RIMAN.asSubstitution(),
                    new Transformer(RenameConflictingIndices.INSTANCE),
                    new Transformer(ExpandBrackets.EXPAND_ALL),
                    IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                    new SqrSubs((SimpleTensor) CC.parse("n_{\\alpha}")),
                    KRONECKER_DIMENSION.asSubstitution(),
                    CalculateNumbers.INSTANCE,
                    new Transformer(CollectPowers.INSTANCE),
                    CollectFactory.createCollectEqualTerms(),
                    CalculateNumbers.INSTANCE
                    ,
                    Averaging.INSTANCE,
                    new Transformer(ExpandBrackets.EXPAND_ALL),
                    IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                    KRONECKER_DIMENSION.asSubstitution(),
                    CalculateNumbers.INSTANCE,
                    new Expression("R_{\\mu \\nu}^{\\mu}_{\\alpha} = R_{\\nu\\alpha}"),
                    IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                    new Expression("R_{\\mu}^{\\mu} = R"),
                    IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                    CollectFactory.createCollectAllEqualTerms(),
                    new Transformer(new POO()),
                    CalculateNumbers.INSTANCE,
                    new Expression("R_{\\mu\\nu}^{\\alpha}_{\\alpha}=0"),
                    //                new Expression("F^{\\alpha}_{\\alpha}=0"),
                    CalculateNumbers.INSTANCE,
                    CollectFactory.createCollectAllEqualTerms(),
                    CalculateNumbers.INSTANCE,
                    new Expression("F_\\mu^\\mu_\\alpha\\beta=0"),
                    new Expression("F_\\alpha\\beta_\\mu^\\mu=0"),
                    IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                    new Expression("R_{\\mu\\nu\\alpha\\beta}*R^{\\mu\\alpha\\nu\\beta}=(1/2)*R_{\\mu\\nu\\alpha\\beta}*R^{\\mu\\nu\\alpha\\beta}"),
                    IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                    new Expression("R_{\\mu\\nu\\alpha\\beta}*R^{\\mu\\nu\\alpha\\beta}=4*R_{\\mu\\nu}*R^{\\mu\\nu}-R*R"),
                    IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                    new Expression("R_{\\mu \\nu}^{\\mu}_{\\alpha} = R_{\\nu\\alpha}"),
                    new Transformer(ExpandBrackets.EXPAND_ALL),
                    IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                    new Expression("R_{\\mu\\nu}^{\\alpha}_{\\alpha}=0"),
                    CollectFactory.createCollectAllEqualTerms(),
                    CalculateNumbers.INSTANCE
                    );
        }

    }

    public final void evalAction() {

        ACTION.eval(
//                Kn,
                //indicesInsertion,
                //                L.asSubstitution(),
                Flat, WR, SR, SSR, FF, FR, RR,
                CalculateNumbers.INSTANCE,
                //RICCI.asSubstitution(),
                //RIMAN.asSubstitution(),
//                new Transformer(RenameConflictingIndices.INSTANCE),
                new Transformer(ExpandBrackets.EXPAND_EXCEPT_SYMBOLS),
                IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                KRONECKER_DIMENSION.asSubstitution());
        ACTION.eval(
//                 Averaging
                new SqrSubs((SimpleTensor) CC.parse("n_{\\alpha}")),
                CalculateNumbers.INSTANCE,
                new Transformer(ExpandBrackets.EXPAND_ALL),
                Averaging.INSTANCE,
                new Transformer(ExpandBrackets.EXPAND_ALL),
                IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                KRONECKER_DIMENSION.asSubstitution(),
                CalculateNumbers.INSTANCE,
                new Expression("R_{\\mu \\nu}^{\\mu}_{\\alpha} = R_{\\nu\\alpha}"),
                IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                new Expression("R_{\\mu}^{\\mu} = R"),
                IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                CollectFactory.createCollectAllEqualTerms(),
                new Transformer(new POO()),
                CalculateNumbers.INSTANCE,
                new Expression("R_{\\mu\\nu}^{\\alpha}_{\\alpha}=0"),
                //                new Expression("F^{\\alpha}_{\\alpha}=0"),
                CalculateNumbers.INSTANCE,
                CollectFactory.createCollectAllEqualTerms(),
                CalculateNumbers.INSTANCE,
//               
                new Expression("F_\\mu\\nu\\alpha\\beta=R_\\mu\\nu\\alpha\\beta"),
                IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                new Expression("R_{\\mu\\nu\\alpha\\beta}*R^{\\mu\\alpha\\nu\\beta}=(1/2)*R_{\\mu\\nu\\alpha\\beta}*R^{\\mu\\nu\\alpha\\beta}"),
                IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                new Expression("R_{\\mu\\nu\\alpha\\beta}*R^{\\mu\\nu\\alpha\\beta}=4*R_{\\mu\\nu}*R^{\\mu\\nu}-R*R"),
                IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                new Expression("R_{\\mu \\nu}^{\\mu}_{\\alpha} = R_{\\nu\\alpha}"),
                new Transformer(ExpandBrackets.EXPAND_ALL),
                IndicesContractionsTransformation.CONTRACTIONS_WITH_METRIC,
                new Expression("R_{\\mu\\nu}^{\\alpha}_{\\alpha}=0"),
                CollectFactory.createCollectAllEqualTerms(),
                CalculateNumbers.INSTANCE);
    }

    private static class POO implements Transformation {
        @Override
        public Tensor transform(Tensor tensor) {
            if (!(tensor instanceof Pow))
                return tensor;
            Pow power = (Pow) tensor;
            if (power.getPower() instanceof TensorNumber && power.getTarget() instanceof TensorNumber) {
                ComplexElement p = ((TensorNumber) power.getPower()).getValue();
                ComplexElement t = ((TensorNumber) power.getTarget()).getValue();
                RationalElement r = t.getReal().pow(p.getReal());
                assert p.getImagine().isZero();
                assert t.getImagine().isZero();
                return new TensorNumber(new ComplexElement(r, RationalElement.ZERO));
            }
            return tensor;
        }
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
            for (Tensor m : MATRIX_ACTION)
                if (TTest.testEqualstensorStructure(tensor, m))
                    return true;
            return false;
        }
    });
}
