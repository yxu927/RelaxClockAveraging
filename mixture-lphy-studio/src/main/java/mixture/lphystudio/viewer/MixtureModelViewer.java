package mixture.lphystudio.viewer;

import lphy.core.model.Value;
import lphystudio.app.graphicalmodelpanel.viewer.Viewer;
import mixture.lphy.evolution.auto.MixturePhyloCTMC;

import javax.swing.*;

public class MixtureModelViewer implements Viewer{
    public MixtureModelViewer(){}

    @Override
    public boolean match(Object value) {
        return value instanceof MixturePhyloCTMC ||
                (value instanceof Value && ((Value) value).value() instanceof MixturePhyloCTMC);

    }

    @Override
    public JComponent getViewer(Object value) {
        if (match(value)) {
            return new JTextArea(value.toString());
        }
        String text = ((Value<MixturePhyloCTMC>) value).value().toString();
        return new JTextArea(text);
    }

    @Override
    public String toString() { return "Read Count Viewer"; }
}

