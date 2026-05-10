package mixture.lphystudio.viewer;

import lphy.core.model.Value;
import lphystudio.app.graphicalmodelpanel.viewer.Viewer;
import mixture.lphy.evolution.auto.MixturePhyloCTMC;

import javax.swing.*;

public class MixtureModelViewer implements Viewer {

    public MixtureModelViewer() {}

    @Override
    public boolean match(Object value) {
        Object obj = unwrap(value);
        return obj instanceof MixturePhyloCTMC;
    }

    @Override
    public JComponent getViewer(Object value) {
        Object obj = unwrap(value);
        if (obj instanceof MixturePhyloCTMC) {
            return new JTextArea(obj.toString());
        }
        return new JTextArea("");
    }

    private Object unwrap(Object value) {
        if (value instanceof Value<?>) {
            return ((Value<?>) value).value();
        }
        return value;
    }

    @Override
    public String toString() {
        return "Mixture Model Viewer";
    }
}