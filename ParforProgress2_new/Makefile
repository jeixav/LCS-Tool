JFLAGS = -g
JC = javac
.SUFFIXES: .java .class
.java.class:
	$(JC) $(JFLAGS) $*.java

CLASSES = \
	ParforProgressServer2.java \
	ParforProgressClient2.java

default: classes

classes: $(CLASSES:.java=.class)

clean:
	$(RM) *.class