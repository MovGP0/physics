package cc.redberry.physics.feyncalc.context;

/**
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public final class FeynCalcContextManager {
    private static FeynCalcContext currentContext;

    public static void setCurrentContext(FeynCalcContext context) {
        currentContext = context;
    }

    public static void setContextToDefault() {
        currentContext = new FeynCalcContext(DefaultFeynCalcSettings.create());
    }

    public static FeynCalcContext getCurrentContext() {
        if (currentContext == null)
            setContextToDefault();
        return currentContext;
    }

}
