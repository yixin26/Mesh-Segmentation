#include "CurveNetMaker.h"
#include <QtWidgets/QApplication>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	CurveNetMaker w;
	w.show();
	return a.exec();
}
