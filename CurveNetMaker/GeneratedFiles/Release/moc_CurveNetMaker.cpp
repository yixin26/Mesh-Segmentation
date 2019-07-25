/****************************************************************************
** Meta object code from reading C++ file 'CurveNetMaker.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.12.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../CurveNetMaker.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'CurveNetMaker.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.12.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_CurveNetMaker_t {
    QByteArrayData data[7];
    char stringdata0[87];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_CurveNetMaker_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_CurveNetMaker_t qt_meta_stringdata_CurveNetMaker = {
    {
QT_MOC_LITERAL(0, 0, 13), // "CurveNetMaker"
QT_MOC_LITERAL(1, 14, 12), // "OpenMeshFile"
QT_MOC_LITERAL(2, 27, 0), // ""
QT_MOC_LITERAL(3, 28, 13), // "SaveCurveFile"
QT_MOC_LITERAL(4, 42, 12), // "AlphaChanged"
QT_MOC_LITERAL(5, 55, 16), // "ShowFeatureLines"
QT_MOC_LITERAL(6, 72, 14) // "SmoothBoundary"

    },
    "CurveNetMaker\0OpenMeshFile\0\0SaveCurveFile\0"
    "AlphaChanged\0ShowFeatureLines\0"
    "SmoothBoundary"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_CurveNetMaker[] = {

 // content:
       8,       // revision
       0,       // classname
       0,    0, // classinfo
       5,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,   39,    2, 0x0a /* Public */,
       3,    0,   40,    2, 0x0a /* Public */,
       4,    0,   41,    2, 0x0a /* Public */,
       5,    0,   42,    2, 0x0a /* Public */,
       6,    0,   43,    2, 0x0a /* Public */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void CurveNetMaker::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        auto *_t = static_cast<CurveNetMaker *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->OpenMeshFile(); break;
        case 1: _t->SaveCurveFile(); break;
        case 2: _t->AlphaChanged(); break;
        case 3: _t->ShowFeatureLines(); break;
        case 4: _t->SmoothBoundary(); break;
        default: ;
        }
    }
    Q_UNUSED(_a);
}

QT_INIT_METAOBJECT const QMetaObject CurveNetMaker::staticMetaObject = { {
    &QMainWindow::staticMetaObject,
    qt_meta_stringdata_CurveNetMaker.data,
    qt_meta_data_CurveNetMaker,
    qt_static_metacall,
    nullptr,
    nullptr
} };


const QMetaObject *CurveNetMaker::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *CurveNetMaker::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_CurveNetMaker.stringdata0))
        return static_cast<void*>(this);
    return QMainWindow::qt_metacast(_clname);
}

int CurveNetMaker::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 5)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 5;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 5)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 5;
    }
    return _id;
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
