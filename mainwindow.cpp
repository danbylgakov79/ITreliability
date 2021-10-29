#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <iomanip>
#include <math.h>

using namespace std;

QVector<QVector<double>> X;
QVector<QVector<double>> Y;
vector<vector<int> > graph;
//set<list<int> > routs;
vector<vector<int>> routs;// все найденные пути между заданной парой вершин
vector<int> path; // найденный путь
set<set<int>> ribs;
int N;
bool used[MAXV];
bool isfile = false;

struct Simulation{
    //время симуляции
    double TimeSimulation;
    //кол-во итераций
    int countIterations;
    //среднееквадратичное для восстановления
    double sigmaRecovery;
    //мат. ожидание для восстановления
    double muRecovery;
    //среднееквадратичное для отказа
    double sigmaFailure;
    //мат. ожидание для отказа
    double muFailure;
};

Simulation sim;

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow) {
    ui->setupUi(this);
    ui->plainTextEdit->setReadOnly(true);
    ui->plainTextEdit_2->setReadOnly(true);
    ui->plainTextEdit_3->setReadOnly(true);
    ui->lineEdit_2->setValidator( new QIntValidator(0, 1000000000, this) );
    ui->lineEdit_3->setValidator( new QIntValidator(0, 1000000000, this) );
    ui->lineEdit_4->setValidator( new QIntValidator(0, 1000000000, this) );
    ui->lineEdit_5->setValidator( new QIntValidator(0, 1000000000, this) );
    ui->lineEdit_6->setValidator( new QIntValidator(0, 1000000000, this) );
    ui->lineEdit_9->setValidator( new QIntValidator(0, 1000000000, this) );
    ui->lineEdit_2->setText(QString::number(876));// мат ожидание для восстановления
    ui->lineEdit_3->setText(QString::number(175)); // среднеквадратичное для восстановления
    ui->lineEdit_4->setText(QString::number(3504));// мат ожидание отказ
    ui->lineEdit_5->setText(QString::number(876));// средне квадрат отказ
    ui->lineEdit_6->setText(QString::number(1000));//итерации
    ui->lineEdit_9->setText(QString::number(43800));//время работы

    sim.muFailure = 3504;
    sim.sigmaFailure = 876;
    sim.muRecovery = 876;
    sim.sigmaRecovery = 175;
    sim.countIterations = 1000;
    sim.TimeSimulation = 43800;

    connect(ui->pushButton,SIGNAL(clicked()),this, SLOT(on_pushButtonDistribution_clicked()));
}

MainWindow::~MainWindow() {
    delete ui;
}

void printVD(string a, QVector<double> arg) {
    if (arg.empty()) {
        cout << "Вектор пуст " << endl;
    }
    else {
        cout << a << endl;
        for (auto it = arg.begin(); it != arg.end(); it++) {
            cout << *it << " ";
        }
        cout << endl;
    }
}

void printVI(string a, QVector<int> arg) {
    if (arg.empty()) {
        cout << "Вектор пуст " << endl;
    }
    else {
        cout << a << endl;
        for (auto it = arg.begin(); it != arg.end(); it++) {
            cout << *it << " ";
        }
        cout << endl;
    }
}

void printVIN(string a, vector<int> arg) {
    if (arg.empty()) {
        cout << "Вектор пуст " << endl;
    }
    else {
        cout << a << endl;
        for (auto it = arg.begin(); it != arg.end(); it++) {
            cout << *it << " ";
        }
        cout << endl;
    }
}

void dfs(int const& nodeCur, int const& nodeLast) {
    path.emplace_back(nodeCur);
    if (nodeCur == nodeLast) {
        routs.push_back(path);
        path.pop_back();
        return;
    }
    for (auto const& val : graph.at(nodeCur - 1)) {
        std::set<int> rib = { nodeCur, val };
        if (ribs.find(rib) == ribs.end()) {
            ribs.emplace(rib);
            dfs(val, nodeLast);
            ribs.erase(ribs.find(rib));
        }
    }
    path.pop_back();
}
void spisSmej(QVector<QVector<int>> g) {
    QVector<QVector<int>> matrixSmej = g;
    for (int i = 0; i < g.size(); i++) {
        for (int j = 0; j < g.size(); j++) {
            if (g[i][j] != 0) {
                matrixSmej[i][j] = 1;
            }
            else {
                matrixSmej[i][j] = 0;
            }
        }
    }
    int n = g.size();
    graph.resize(n + 1);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            if (matrixSmej[i][j] != 0) graph[i].push_back(j + 1);
        }
}

void searchWays(QVector<QVector<double>> MatrixElements) {
    int x = (int)MatrixElements[0][0];
    int y = 0;
    for(int i = 0; i < MatrixElements.size(); i++) {
        for(int j = 0;j < MatrixElements[i].size(); j++) {
            y++;
        }
    }
    pair<int, int> nodes(x, y);
    dfs(nodes.first, nodes.second);
}

void checkWays(QVector<int>& ways, vector<vector<int>> routs, int x) {
    int min = 100000;
    QVector<int> waysTemp;
    for (int i = 0; i < routs.size(); i++) {
        if (routs[i].size() < min) {
            min = routs[i].size();
        }
    }
    for (int i = 0; i < routs.size() ; i++) {
        if (routs[i].size() == min) {
            waysTemp.push_back(i);
        }
    }

   //int k=0;
   int sum=0;
  ways.push_back(waysTemp[0]);
  // while(k!=waysTemp.size()-1)
  for(int k=1;k<waysTemp.size();k++)
   {
       //int wk = waysTemp[k];
       for(int i=k;i<waysTemp.size();i++)
       {


           for(int j=0;j<routs[waysTemp[i]].size();j++)
           {
              if (routs[waysTemp[k]][j]!=routs[waysTemp[i]][j])
              {
                  sum++;
              }
           }
           if(sum>=x)
           {
               k=i;

               ways.push_back(waysTemp[i]);
               sum=0;

           }
       }
       sum=0;
   }


}

void printVVD(string a, QVector<QVector<double>> arg) {
    if (arg.empty()) {
        cout << "Вектор пуст " << endl;
    }
    else {
        cout << a << "[" << arg.size() << "] = " << endl;
        for (auto it = arg.begin(); it != arg.end(); ++it) {
            for (const auto& i : *it)
                cout << i << " ";
            cout << endl;
        }
        cout << endl;
    }
}

void printVVI(string a, vector<vector<int>> arg) {
    if (arg.empty()) {
        cout << "Вектор пуст " << endl;
    }
    else {
        cout << a << "[" << arg.size() << "] = " << endl;
        for (auto it = arg.begin(); it != arg.end(); ++it) {
            for (const auto& i : *it)
                cout << i << " ";
            cout << endl;
        }
        cout << endl;
    }
}

void removevecV(int countCoeffOfElement, QVector<double>& a) {
    QVector<double> temp;
    if (countCoeffOfElement == 0) return;
    for (int i = 0; i < countCoeffOfElement; i++) {
        temp.push_back(a[i]);
    }
    a.clear();
    QVector<double>().swap(a);
    for (int i = 0; i < temp.size(); i++) {
        a.push_back(temp[i]);
    }
}

void removevecVV(int countCoeffOfElement, QVector<QVector<double>>& a) {
    QVector<QVector<double>> temp;
    if (countCoeffOfElement == 0) return;
    for (int i = 0; i < countCoeffOfElement; i++) {
        temp.push_back(a[i]);
    }
    a.clear();
    QVector<QVector<double>>().swap(a);
    for (int i = 0; i < temp.size(); i++) {
        a.push_back(temp[i]);
    }
}

//Генерация случайных чисел полярным методом
double NormalDistribution(double mean, double stdDev) {
    static double spare;
    static bool hasSpare = false;
    if (hasSpare) {
        hasSpare = false;
        return spare * stdDev + mean;
    }
    else {
        double u, v, s;
        do {
            u = (rand() / ((double)RAND_MAX)) * 2.0 - 1.0;
            v = (rand() / ((double)RAND_MAX)) * 2.0 - 1.0;
            s = u * u + v * v;
        } while (s >= 1.0 || s == 0.0);
        s = sqrt(-2.0 * log(s) / s);
        spare = v * s;
        hasSpare = true;
        return mean + stdDev * u * s;
    }
}

void CreateGraphicDistribution(QVector<double> &Xx,QVector<double> &Yy) {
    int N = 100;
    double dX = (2 * sim.muFailure) / N;
    for (int j = 0; j < N; j++) {
        Xx.push_back(j * dX);
    }
    for (int i = 0; i < Xx.size(); i++) {
        Yy.push_back((1 / (sim.sigmaFailure * sqrt(2 * acos(-1)))) * exp((-pow(Xx[i] - sim.muFailure, 2)) / pow(sim.sigmaFailure, 2)));
    }
}

void readFile(QString fname, QVector<QVector<int>>& W) {
    std::string filename = fname.toUtf8().constData();
    filename = filename + ".txt";
    ifstream file(filename);
    if(!file) {
        isfile = true;
        return;
    }
    int N = 0;
    string str;
    //считаем кол-во строк
    while (getline(file, str)) {
        N++;
    }
    int M = 0;
    file.clear();
    file.seekg(0);
    int t = 0;
    while (file >> t) {
        M++;
    }
    M /= N;
    //возвращаемся к началу файла
    file.clear();
    file.seekg(0);
    W.resize(N);
    for (int i = 0; i < W.size(); i++) {
        W[i].resize(M);
        for (int j = 0; j < W[i].size(); j++) {
            file >> W[i][j];
        }
    }
}

void MainWindow::on_pushButtonDistribution_clicked() {
    ui->plainTextEdit->clear();
    QVector<double> Xx;
    QVector<double> Yy;
    CreateGraphicDistribution(Xx,Yy);
    ui->plainTextEdit->setPlainText("Параметры графика нормального распределения\n");
    ui->plainTextEdit->setPlainText(ui->plainTextEdit->toPlainText() + "Мат. ожидание = " + QString::number(sim.muRecovery) + "\n");
    ui->plainTextEdit->setPlainText(ui->plainTextEdit->toPlainText() + "Сред. квадр. отклонение = " + QString::number(sim.muFailure) + "\n");
    double a = Xx[0]; //Начало интервала, где рисуем график по оси Ox
    double b = Xx[Xx.size() - 1]; //Конец интервала, где рисуем график по оси Ox
    //double h = 0.01; //Шаг, с которым будем пробегать по оси Ox
    int N = Xx.size(); //Вычисляем количество точек, которые будем отрисовывать
    ui->widget->clearGraphs();//Если нужно, но очищаем все графики
    //Добавляем один график в widget
    ui->widget->addGraph();
    //Говорим, что отрисовать нужно график по нашим двум массивам x и y
    ui->widget->graph(0)->setData(Xx, Yy);
    //Подписываем оси Ox и Oy
    ui->widget->xAxis->setLabel("x");
    ui->widget->yAxis->setLabel("y");
    //Установим область, которая будет показываться на графике
    ui->widget->xAxis->setRange(a, b);//Для оси Ox
    //Для показа границ по оси Oy сложнее, так как надо по правильному
    //вычислить минимальное и максимальное значение в векторах
    double minY = Yy[0], maxY = Yy[0];
    for (int i = 1; i < N; i++) {
        if (Yy[i] < minY) minY = Yy[i];
        if (Yy[i] > maxY) maxY = Yy[i];
    }
    ui->widget->yAxis->setRange(minY, maxY);//Для оси Oy
    //И перерисуем график на нашем widget
    ui->widget->replot();
}

void time(QVector<QVector<double>>& States) {
    double ItemIsIdle = 0;//Элемент в восстановлении
    double Itemserves = 0;//Элемент в отказе
    double valueOfWorkingStates = 0; // значение длительности времени восстановления
    double valueOfNonWorkingStates = 0; // значение длительности времени отказа
    States.push_back(QVector<double>());
    // 1 - восстановление
    // 0 - отказ
    //сначала элемент идет в отказ, то есть идет состояние восстановления
    States[0].push_back(1);
    //получаем значение
    Itemserves = NormalDistribution(sim.muFailure, sim.sigmaFailure);
    ItemIsIdle = Itemserves;
    //записываем время работы до отказа
    int countOfStates = 0;
    double mainTime = 0;
    while(mainTime < sim.TimeSimulation) {
        mainTime += 1;
        if (States[countOfStates][0] == 0) {
            countOfStates++;
            States.push_back(QVector<double>());
            States[countOfStates].push_back(1);
        }
        if (mainTime >= ItemIsIdle) {
            valueOfWorkingStates = NormalDistribution(sim.muRecovery, sim.sigmaRecovery);
            valueOfNonWorkingStates = Itemserves;
            //длительность времени работы до отказа
            States[countOfStates].push_back(valueOfNonWorkingStates);
            countOfStates++;
            States.push_back(QVector<double>());
            States[countOfStates].push_back(0);
            //длительность времени восстановления
            States[countOfStates].push_back(valueOfWorkingStates);
            //время работы до отказа
            Itemserves = NormalDistribution(sim.muFailure, sim.sigmaFailure);
            //время работы до отказа
            ItemIsIdle += Itemserves;
            mainTime += valueOfWorkingStates;
        }
    }
    States[countOfStates].push_back(Itemserves);
}


double FunctNormal(double x,double mu,double sigma) {
    double y;
    //double x = NormalDistribution(mu,sigma);
   y = ((1 / (sigma * sqrt(2 * acos(-1)))) * exp((-pow((x - mu), 2)) / pow(sigma, 2)));
    //y= (x*x)/2;
   // cout<<"y="<<y;
    return y;
}

const double eps = 0.1e-5;
double func(double x) {
    return std::exp(-x * x / 2);
}

double simpsonIntegr(double a, double b) {
    return ((b - a) / 6) * (func(a) + 4 * func((a + b) / 2) + func(b));
}
double FUN(double x)
{
 return exp(-x * x / 2);
}

double QUANC8(double A, double B, double ABSERR, double RELERR,
            double RESULT, double ERREST, int NOFUN, double FLAG)
{ double W0, W1, W2, W3, W4, AREA, X0, F0, STONE, STEP, COR11, TEMP;
  double QPREV, QNOW, QDIFF, QLEFT, ESTERR, TOLERR;
  double QRIGHT[32], F[17], X[17], FSAVE[9][31], XSAVE[9][31];
  int LEVMIN, LEVMAX, LEVOUT, NOMAX, NOFIN, LEV, NIM, I, J;

  RELERR = 1.0E-10;
  ABSERR = 0.0;

// *** 芬일 1 ***
  LEVMIN = 1;
  LEVMAX = 30;
  LEVOUT = 6;
  NOMAX = 5000;
  NOFIN = NOMAX-8*(LEVMAX - LEVOUT + pow(2,LEVOUT + 1));
  W0 = 3956.0/14175.0;
  W1 = 23552.0/14175.0;
  W2 = -3712.0/14175.0;
  W3 = 41984.0/14175.0;
  W4 = -18160.0/14175.0;
  FLAG = 0.0;
  RESULT = 0.0;
  COR11 = 0.0;
  ERREST = 0.0;
  AREA = 0.0;
  NOFUN = 0;
  if (A == B) return RESULT;

// ***芬일 2***
  LEV = 0;
  NIM = 1;
  X0 = A;
  X[16] = B;
  QPREV = 0.0;
  F0 = FUN(X0);
  STONE = (B - A)/16.0;
  X[8] = (X0 + X[16])/2.0;
  X[4] = (X0 + X[8])/2.0;
  X[12] = (X[8]+X[16])/2.0;
  X[2] = (X0+X[4])/2.0;
  X[6] = (X[4]+X[8])/2.0;
  X[10] = (X[8]+X[12])/2.0;
  X[14] = (X[12]+X[16])/2.0;
  for (J=2; J<=16; J=J+2)
     F[J]=FUN(X[J]);
  NOFUN = 9;

// ***芬일 3***
m30: X[1] = (X0 + X[2])/2.0;
  F[1] = FUN(X[1]);
  for (J=3; J<=15; J=J+2)
    { X[J] = (X[J-1] + X[J+1])/2.0;
      F[J] = FUN(X[J]);
    }
  NOFUN = NOFUN + 8;
  STEP = (X[16] - X0)/16.0;
  QLEFT = (W0*(F0+F[8])+W1*(F[1]+F[7])+W2*(F[2]+F[6])+W3*(F[3]+F[5])+W4*F[4])*STEP;
  QRIGHT[LEV+1] = (W0*(F[8]+F[16])+W1*(F[9]+F[15])+W2*(F[10]+F[14])+W3*(F[11]+F[13])+W4*F[12])*STEP;
  QNOW=QLEFT+QRIGHT[LEV+1];
  QDIFF = QNOW-QPREV;
  AREA = AREA+QDIFF;

// ***芬일 4***
  ESTERR = fabs(QDIFF)/1023.0;
  TOLERR = max(ABSERR, RELERR*fabs(AREA))*(STEP/STONE);
  if (LEV < LEVMIN) goto m50;
  if (LEV >= LEVMAX) goto m62;
  if (NOFUN > NOFIN) goto m60;
  if (ESTERR <= TOLERR) goto m70;

// ***芬일 5***
m50: NIM = 2*NIM;
  LEV = LEV + 1;
  for (I=1; I<=8; I++)
    { FSAVE[I][LEV] = F[I+8];
      XSAVE[I][LEV] = X[I+8];
    }
  QPREV = QLEFT;
  for (I=1; I<=8; I++)
    { J = - I;
      F[2*J+18] = F[J+9];
      X[2*J+18] = X[J+9];
    }
  goto m30;

// ***芬일 6***
m60: NOFIN = 2*NOFIN;
  LEVMAX = LEVOUT;
  FLAG=FLAG+(B - X0)/(B-A);
  goto m70;
m62: FLAG = FLAG + 1.0;

// ***芬일 7***
m70: RESULT = RESULT + QNOW;
  ERREST = ERREST + ESTERR;
  COR11 = COR11 + QDIFF/1023.0;
m72:
  if (NIM == 2*(NIM/2)) goto m75;
  NIM = NIM/2;
  LEV = LEV - 1;
  goto m72;
m75: NIM = NIM + 1;
  if (LEV <= 0) goto m80;
  QPREV = QRIGHT[LEV];
  X0 = X[16];
  F0 = F[16];
  for (I=1; I<=8; I++)
    { F[2*I] = FSAVE[I][LEV];
      X[2*I] = XSAVE[I][LEV];
    }
  goto m30;

// ***芬일 8***
m80: RESULT = RESULT + COR11;
  if (ERREST == 0.0) return RESULT;
m82: TEMP = fabs(RESULT) + ERREST;
  if (TEMP != fabs(RESULT)) return RESULT;
  ERREST = 2.0*(ERREST);
  goto m82;
}

double LaplasFunction(double x)
{ double A, B, ABSERR, RELERR, RESULT, ERREST, FLAG;
  int NOFUN;

  double a = 0.0;
  double b = x;
  double tst, tst1, incr;
  unsigned int c = 2;
//  tst1 = simpsonIntegr(a, b);
  A=a; B=b;
  tst1 = QUANC8(A, B, ABSERR, RELERR, RESULT, ERREST, NOFUN, FLAG);
  do { tst = tst1;
       tst1 = 0;
       incr = (b - a) / c;
       for (unsigned int i = 0; i < c; ++i)
         { //        tst1 = tst1 + simpsonIntegr(a + incr * i, a + incr * (i + 1));
//           A=a; B=b;
           tst1 = tst1 + QUANC8(A+incr*i, A+incr*(i+1), ABSERR, RELERR, RESULT, ERREST, NOFUN, FLAG);
         }
       c += 1;
     } while (fabs(tst - tst1) >= eps);
  tst1 *= (1. / sqrt(2. * M_PI));
  return tst1;
}/*
double LaplasFunction(double x, double eps = ::eps) {
    double a = 0.0;
    double b = x;
    double tst, tst1, incr;
    unsigned int c = 2;
    tst1 = simpsonIntegr(a, b);
    do {
        tst = tst1;
        tst1 = 0;
        incr = (b - a) / c;
        for (unsigned int i = 0; i < c; ++i) {
            tst1 = tst1 + simpsonIntegr(a + incr * i, a + incr * (i + 1));
        }
        c += 1;
    } while (fabs(tst - tst1) >= eps);
    tst1 *= (1. / sqrt(2. * M_PI));
    return tst1;
}*/
void CoeffAnalytic(QVector<double>& Coeff,double num){
    double t=sim.TimeSimulation;
        double T=sim.muFailure, Tb=sim.muRecovery,o=sim.sigmaFailure,ob=sim.sigmaRecovery;
        double number = num;
        double s = t/number;
        //float x;
        Coeff.push_back(1);

        float start=0+s;
        for(double d =start;d<t;d+=s){
            double temp=0;
        for(double i=1;i<number;i++)
        {
             temp += LaplasFunction((d-(i*T)-(i*Tb))/sqrt((i*pow(o,2))+(i*pow(ob,2))))-LaplasFunction((d-((i+1)*T)-(i*Tb))/sqrt(((i+1)*pow(o,2))+(i*pow(ob,2))));
           // cout<<"tm ="<<temp<<endl;
        }
           // cout<<"tm3 ="<<temp<<endl;
            double tm = (1./2.)-LaplasFunction((d-T)/o);
            //cout<<"tm2 ="<<tm<<endl;
            Coeff.push_back(tm+temp);


        }
}

//функция, рассчитывающая коэффициент готовности для одного элемента
double calculateCoeffOfElement(QVector<double>& CoeffsOfElement, QVector<double>& Recovery, QVector<double>& Failure, QVector<double>& BinaryElement) {
    double valueOfWorkingStates = 0; //сумма значений рабочих состояний
    double valueOfNonWorkingStates = 0; //сумма значений не рабочих состояний
    double sumOfStates; //общая сумма состояний
    QVector<QVector<double>> States; //вектор, хранящий состояния (1 - восстановление, 0 - отказ)
    int x = 250;
    for (int i = 0; i < x; i++) {
        if(i%2 == 0){
            States.push_back(QVector<double>());
            States[i].push_back(1);
            States[i].push_back(0);
        }
        else {
            States.push_back(QVector<double>());
            States[i].push_back(0);
            States[i].push_back(0);
        }
    }
    for (int i = 0; i < x; i++) {
        CoeffsOfElement.push_back(0);
    }
    int mincount = x;
    QVector<QVector<double>> TempStates;
    QVector<double> TempCoeffOfElement;

    for (int i = 0; i < sim.countIterations; i++) {
        QVector<double>().swap(TempCoeffOfElement);
        QVector<QVector<double>>().swap(TempStates);
        valueOfWorkingStates = 0;
        valueOfNonWorkingStates = 0;
        time(TempStates);
        for (int j = 0; j < TempStates.size();) {
            if (TempStates[j][0] != 0 && j < TempStates.size()) {
                valueOfWorkingStates += TempStates[j][1];
                j++;
            }
            /*Иначе если текущее состояние = 1, то суммируем время восстановиления*/
            else if (TempStates[j][0] != 1 && j < TempStates.size()) {
                valueOfNonWorkingStates += TempStates[j][1];
                j++;
            }
            sumOfStates = valueOfNonWorkingStates + valueOfWorkingStates;
            TempCoeffOfElement.push_back(valueOfWorkingStates / sumOfStates);
        }
        if (mincount > TempStates.size()){
            mincount = TempStates.size();
        }
        for (int k = 0; k < mincount; k++) {
            States[k][1] += TempStates[k][1];
        }
        for (int k = 0; k < mincount; k++) {
            CoeffsOfElement[k] += TempCoeffOfElement[k];
        }
    }

    for (int i = 0; i < mincount; i++) {
        States[i][1] /= sim.countIterations;
    }
    for (int i = 0; i < mincount; i++) {
        CoeffsOfElement[i] /= sim.countIterations;
    }
    removevecV(mincount, CoeffsOfElement);
    removevecVV(mincount, States);

    for (int i = 0; i < States.size(); i++) {
        if (i % 2 == 0) {
            Failure.push_back(States[i][1]);
        }
        else {
            Recovery.push_back(States[i][1]);
        }
    }

    double sum = 0;
    for(int i = 0; i < States.size(); i++) {
        sum += States[i][1];
    }

    for (int i = 0; i < States.size(); i++) {
            if (States[i][0] == 1) {
                for (int j = 0; j < (int)(States[i][1] * 60); j++) {
                    BinaryElement.push_back(1);
                }
            }
            else {
                for (int j = 0; j < (int)(States[i][1] * 60); j++) {
                    BinaryElement.push_back(0);
                }
            }
    }
    sum = 0;


    double TotalCoeffOfElement = CoeffsOfElement.back();
    return TotalCoeffOfElement;
}

int minimalSizeBinary(QVector<QVector<double>> BinaryElements) {
    int minSize = BinaryElements[0].size();
    for(int i = 1; i < BinaryElements.size(); i++) {
        if(BinaryElements[i].size() < minSize) {
            minSize = BinaryElements[i].size();
        }
    }
    return minSize;
}

//Расчет иммитационного коэффа верхняя граница
double TotalCoeffSystemImitationUpper(QVector<QVector<double>> &BinaryElements, QVector<int> ways,QVector<QVector<double>> W) {
    int minSize = minimalSizeBinary(BinaryElements);
    QVector<QVector<double>> BinaryForCoeff(routs[ways[0]].size());
    QVector<double> temp(minSize);

    for (int m = 0; m < minSize; m++) {
        temp[m] = BinaryElements[0][m];
    }

   /* for (int k = 0; k < routs[ways[0]].size(); k++) {

        for (int i = 0; i <= ways.size(); i++) {
            cout<<routs[ways[i]][k] - 1;
            for (int j = 0; j < minSize; j++) {
                temp[j] = temp[j] || BinaryElements[routs[ways[i]][k] - 1][j];
            }
        }
        BinaryForCoeff[k] = temp;

    }
    for (int m = 0; m < minSize; m++) {
        temp[m] = BinaryElements[0][m];
    }
*/

    for (int k = 0; k < routs[ways[0]].size(); k++) {


        for (int i = 0; i < ways.size(); i++) {


            for (int j = 0; j < minSize; j++) {
                temp[j] = temp[j] || BinaryElements[routs[ways[i]][k] - 1][j];
            }
        }

        BinaryForCoeff[k] = temp;

    }

    QVector<double> TotalBinary(minSize);
    QVector<double> temp2 = BinaryForCoeff[0];

    for(int i = 1; i < BinaryForCoeff.size(); i++) {
        for(int j = 0; j < BinaryForCoeff[i].size(); j++) {
            temp2[j] = temp2[j] && BinaryForCoeff[i][j];
        }
    }

    TotalBinary = temp2;
    double countOfWork = 0;
    for (int i = 0; i < TotalBinary.size(); i++) {
        if(TotalBinary[i] == 1) countOfWork++;
    }

    double tac = countOfWork/TotalBinary.size();

    //double totalCoeff = 1 - pow((1 - pow(tac, W.size())), W[0].size()) ;
    double totalCoeff = pow((1 - pow((1.- tac),W.size())),W[0].size());
    double step = TotalBinary.size()/100;
    double k = 0;

    double sum = 0;
    for(int i = 0; i <TotalBinary.size() - step; i+=step) {
        for(int j = i; j <= (i + step); j++) {
            if(TotalBinary[j] == 1) sum++;
        }
        X[3].push_back(k);
        double x = sum/TotalBinary.size();
        double m = pow((1 - pow((1.- x),W.size())),W[0].size());
        Y[3].push_back(m);
        m = 0;
        k += step;
    }
    return totalCoeff;
}
//Расчет нижней границы коэффициента готовности
double TotalCoeffSystemImitationLower(QVector<QVector<double>> &BinaryElements,QVector<int> ways,QVector<double>& KrForWay, QVector<QVector<double>> W) {
    int minSize = minimalSizeBinary(BinaryElements);
    QVector<QVector<double>> BinaryForCoeff(ways.size());
    QVector<double> temp(minSize);

    for (int m = 0; m < minSize; m++) {
        temp[m] = BinaryElements[0][m];
    }

    for (int k = 0; k < ways.size(); k++) {


        for (int i = 0; i < routs[ways[k]].size(); i++) {


            for (int j = 0; j < minSize; j++) {
                temp[j] = temp[j] && BinaryElements[routs[ways[k]][i] - 1][j];
            }
        }

        BinaryForCoeff[k] = temp;

    }

    double sumKr = 0;
    for(int i = 0; i < BinaryForCoeff.size(); i++) {
        for(int j = 0; j < BinaryForCoeff[i].size(); j++) {
            if(BinaryForCoeff[i][j] == 1) sumKr++;
        }
        KrForWay.push_back(sumKr/BinaryForCoeff[i].size());
        sumKr = 0;
    }

    QVector<double> TotalBinary(minSize);
    QVector<double> temp2 = BinaryForCoeff[0];

    for(int i = 1; i < BinaryForCoeff.size(); i++) {
        for(int j = 0; j < BinaryForCoeff[i].size(); j++) {
            temp2[j] = temp2[j] || BinaryForCoeff[i][j];
        }
    }

    TotalBinary = temp2;
    double countOfWork = 0;
    for (int i = 0; i < TotalBinary.size(); i++) {
        if(TotalBinary[i] == 1) countOfWork++;
    }

    double tac = countOfWork/TotalBinary.size();


   // double totalCoeff = pow((1 - pow((1.- tac),W.size())),W[0].size());
    double totalCoeff = 1 - pow((1 - pow(tac, W.size())), W[0].size()) ;
    double step = TotalBinary.size()/100;
    double k = 0;
    double sum = 0;
    for(int i = 0; i < TotalBinary.size() - step; i+=step) {
        for(int j = i; j <= (i + step); j++) {
            if(TotalBinary[j] == 1) sum++;
        }
        X[4].push_back(i);
        double x = sum/TotalBinary.size();
        double m = 1 - pow((1 - pow(x, W.size())), W[0].size()) ;
        Y[4].push_back(m);
        m = 0;
        k++;
    }
    cout<<"check*";
    return totalCoeff;
}

void interpol1(QVector<double> &x, QVector<double> &y) {
    int n = X[0].size();

    QVector<double> a(n);
    QVector<double> b(n);
    QVector<double> c(n);
    QVector<double> h(n);
    QVector<double> d(n);
    QVector<double> k(n);
    QVector<double> l(n);
    QVector<double> r(n);
    QVector<double> s(n);
    //============================нахождение коэффициентов============//

    for (int i = 2; i <= n - 1; i++) {
        k[1] = 0;
        l[1] = 0;

        h[i - 1] = X[0][i - 1] - X[0][i - 2];
        h[i] = X[0][i] - X[0][i - 1];

        s[i] = 2 * (h[i] + h[i - 1]);

        r[i] = 3 * ((Y[0][i] - Y[0][i - 1]) / h[i] - (Y[0][i - 1] - Y[0][i - 2]) / h[i - 1]);

        k[i] = (r[i] - h[i - 1] * k[i - 1]) / (s[i] - h[i - 1] * l[i - 1]);

        l[i] = h[i] / (s[i] - h[i - 1] * l[i - 1]);
    }


    c[n - 1] = k[n - 1];
    for (int i = n - 2; i >= 2; i--)
        c[i] = k[i] - l[i] * c[i + 1];


    for (int i = 1; i <= n - 2; i++) {
        h[i] = X[0][i] - X[0][i - 1];
        a[i] = Y[0][i - 1];
        b[i] = (Y[0][i] - Y[0][i - 1]) / h[i] - h[i] * (2 * c[i] + c[i + 1]) / 3.;
        d[i] = (c[i + 1] - c[i]) / (3. * h[i]);
    }
    int i = 1;
    double x1 = X[0][0];
    double x2;
    double y2 = 0;
    do {
        do {

            y2 = a[i] + b[i] * (x1 - X[0][i - 1]) + c[i] * pow((x1 - X[0][i - 1]), 2) + d[i] * pow((x1 - X[0][i - 1]), 3);
            x.push_back(x1);
            y.push_back(y2);

            x1 = x1 + 1;
            x2 = static_cast<double> (x1);
        } while (x2 < X[0][i]);

        i++;
        x1 = X[0][i - 1];
    } while (i != n - 1);
}

void GrafickCoeff(QVector<double> coeff, double timeSimulation) {
        double g = timeSimulation / (double)coeff.size();
        double j = g;
        int p = 0;
        double b = 0;
        for (int i = 0; i < coeff.size(); i++) {
            b += coeff[i];
        }
        b = b / (coeff.size());
        b -= 0.0001;
        for (int i = coeff.size() - 2; i > 0; i--) {
            if (abs(b - coeff[i]) >= 0.4) {
                goto Leave;
            }
            else {
                if (abs(b - coeff[i]) >= 0.0004) {
                    break;
                }
                else {
                    goto Leave;
                }
            }
        }
        for (int i = coeff.size() - 1; i > 0; i--) {
            if (coeff[i] <= b) {
                coeff.erase(coeff.begin()+i);
            }
            else {
                break;
            }
        }
    Leave:
        if (timeSimulation >= 10000) {

           for (int i = 0; i < coeff.size(); i++) {
                if (i != coeff.size() - 5) {
                   // X[0].push_back(j);
                    Y[0].push_back(coeff[i]);
                    j += g;
                }
                else {
                    while (p != 5) {
                        b = trunc(coeff[i]*1000);
                        b = b / 1000;
                        //X[0].push_back(j);
                        Y[0].push_back(b);
                        j = j + g;
                        p++;
                    }
                    break;
                }
            }
        }
        else {
            for (int i = 0; i < coeff.size(); i++) {
               // X[0].push_back(j);
                Y[0].push_back(coeff[i]);
                j += g;
            }
        }

    double x0=0;
    X[0].push_back(x0);
    x0+=sim.muFailure;
    for(int i=1;i<Y[0].size()-1;i++)
    {
    if(Y[0][i]>Y[0][i+1])
    {
        //x0+=sim.muRecovery;

        X[0].push_back(x0);
        x0+=sim.muFailure;
    }
    else {

        X[0].push_back(x0);
        x0+=sim.muRecovery;

    }
    }
    QVector<double> x1;
    QVector<double> y2;

    interpol1(x1,y2);

    X[0].swap(x1);
    Y[0].swap(y2);
    /*double x=X[0][X[0].size()-1];
    while(Y[0].size()!=Y[6].size()) {
        Y[0].push_back(Y[0][Y[0].size()-1]);
      //  Y[0].push_back(Y[0][Y[0].size()-1]);
        x+=sim.muFailure;
        X[0].push_back(x);
      //  X[0].push_back(sim.muRecovery);

    }*/

}



void interpol(QVector<double> &x, QVector<double> &y) {
    int n = X[6].size();

    QVector<double> a(n);
    QVector<double> b(n);
    QVector<double> c(n);
    QVector<double> h(n);
    QVector<double> d(n);
    QVector<double> k(n);
    QVector<double> l(n);
    QVector<double> r(n);
    QVector<double> s(n);
    //============================нахождение коэффициентов============//

    for (int i = 2; i <= n - 1; i++) {
        k[1] = 0;
        l[1] = 0;

        h[i - 1] = X[6][i - 1] - X[6][i - 2];
        h[i] = X[6][i] - X[6][i - 1];

        s[i] = 2 * (h[i] + h[i - 1]);

        r[i] = 3 * ((Y[6][i] - Y[6][i - 1]) / h[i] - (Y[6][i - 1] - Y[6][i - 2]) / h[i - 1]);

        k[i] = (r[i] - h[i - 1] * k[i - 1]) / (s[i] - h[i - 1] * l[i - 1]);

        l[i] = h[i] / (s[i] - h[i - 1] * l[i - 1]);
    }


    c[n - 1] = k[n - 1];
    for (int i = n - 2; i >= 2; i--)
        c[i] = k[i] - l[i] * c[i + 1];


    for (int i = 1; i <= n - 2; i++) {
        h[i] = X[6][i] - X[6][i - 1];
        a[i] = Y[6][i - 1];
        b[i] = (Y[6][i] - Y[6][i - 1]) / h[i] - h[i] * (2 * c[i] + c[i + 1]) / 3.;
        d[i] = (c[i + 1] - c[i]) / (3. * h[i]);
    }
    int i = 1;
    double x1 = X[6][0];
    double x2;
    double y2 = 0;
    do {
        do {

            y2 = a[i] + b[i] * (x1 - X[6][i - 1]) + c[i] * pow((x1 - X[6][i - 1]), 2) + d[i] * pow((x1 - X[6][i - 1]), 3);
            x.push_back(x1);
            y.push_back(y2);

            x1 = x1 + 1;
            x2 = static_cast<double> (x1);
        } while (x2 < X[6][i]);

        i++;
        x1 = X[6][i - 1];
    } while (i != n - 1);
}


/*void GrafickCoeffA(QVector<double> coeff, double timeSimulation) {
    double g = timeSimulation / (double)coeff.size();
    double j = g;
    int p = 0;
    double b = 0;
    for (int i = 0; i < coeff.size(); i++) {
        b += coeff[i];
    }
    b = b / (coeff.size());
    b -= 0.0001;
    for (int i = coeff.size() - 2; i > 0; i--) {
        if (abs(b - coeff[i]) >= 0.04) {
            goto Leave;
        }
        else {
            if (abs(b - coeff[i]) >= 0.0004) {
                break;
            }
            else {
                goto Leave;
            }
        }
    }
    for (int i = coeff.size() - 1; i > 0; i--) {
        if (coeff[i] <= b) {
            coeff.erase(coeff.begin()+i);
        }

        else {
            break;
        }
    }
Leave:
    if (timeSimulation >= 10000) {

       for (int i = 0; i < coeff.size(); i++) {
            if (i != coeff.size() - 5) {
                X[6].push_back(j);
                Y[6].push_back(coeff[i]);
                j += g;
            }
            else {
                while (p != 5) {
                    b = trunc(coeff[i]*1000);
                    b = b / 1000;
                    X[6].push_back(j);
                    Y[6].push_back(b);
                    j = j + g;
                    p++;
                }
                break;
            }
        }
    }
    else {
        for (int i = 0; i < coeff.size(); i++) {
            X[6].push_back(j);
            Y[6].push_back(coeff[i]);
            j += g;
        }
    }

    QVector<double> x1;
    QVector<double> y2;

    interpol(x1,y2);

    X[6].swap(x1);
    Y[6].swap(y2);
}
*/
void GrafickCoeffA(QVector<double> coeff, double timeSimulation) {
    double g = timeSimulation / (double)coeff.size();
    double j = g;
    int p = 0;
    double b = 0;
    for(int i=0;i<coeff.size();i++)
    {
        if(coeff[i]>1.)
        {
            coeff.erase(coeff.begin()+i);
        }
    }
    for (int i = 0; i < coeff.size(); i++) {
        b += coeff[i];
    }
    b = b / (coeff.size());
    b -= 0.0001;

    for (int i = coeff.size() - 2; i > 0; i--) {
        if (abs(b - coeff[i]) >= 0.04) {
            goto Leave;
        }
        else {
            if (abs(b - coeff[i]) >= 0.001) {
                break;
            }
            else {
                goto Leave;
            }
        }
    }
    for (int i = coeff.size() - 1; i > 0; i--) {
        if (coeff[i] <= b) {
            coeff.erase(coeff.begin()+i);
        }

        else {
            break;
        }
    }
Leave:
    if (timeSimulation >= 10000) {

       for (int i = 0; i < coeff.size(); i++) {
            if (i != coeff.size() - 5) {
                X[6].push_back(j);
                Y[6].push_back(coeff[i]);
                j += g;
            }
            else {
                while (p != 5) {
                    b = trunc(coeff[i]*1000);
                    b = b / 1000;
                    X[6].push_back(j);
                    Y[6].push_back(b);
                    j = j + g;
                    p++;
                }
                break;
            }
        }
    }
    else {
        for (int i = 0; i < coeff.size(); i++) {
            X[6].push_back(j);
            Y[6].push_back(coeff[i]);
            j += g;
        }
    }

    QVector<double> x1;
    QVector<double> y2;

    interpol(x1,y2);

    X[6].swap(x1);
    Y[6].swap(y2);
}

//Вывод в консоль структуры
QString printStruct(QVector<QVector<double>> MatrixElements) {
    string out = "Структура модели  \n";
    int x = 0;
    for(int i = 0; i < MatrixElements.size(); i++) {
        for(int j = 0; j < MatrixElements.size(); j++) {
            x++;
            out += to_string(x)+" ";
        }
        out += "\n";
    }
    QString ent = QString::fromStdString(out);
    return ent;
}
//Заполнение вектора графика времени
void GrafickTime(QVector<QVector<double>>& resultsRecovery,QVector<QVector<double>>& resultsFailure, QVector<QVector<double>> coeff,int k) {
    double x = 0.;
    double timeSum = 0;
    for(int j = 0; j < resultsRecovery[k].size(); j++) {
        timeSum +=resultsRecovery[k][j];
    }
    for(int j = 0; j < resultsFailure[k].size(); j++) {
        timeSum +=resultsFailure[k][j];
       }
    Y[1].push_back(coeff[k][0]);
    Y[1].push_back(coeff[k][1]);
    X[1].push_back(x);
    x += resultsFailure[k][0];
    X[1].push_back(x);
    int l = 0;
    int m = 1;
    for(int i = 2 ;i < coeff[k].size(); i++) {
        if(i%2 == 0) {
           Y[1].push_back(coeff[k][i]);
           Y[1].push_back(coeff[k][i]);
           X[1].push_back(x);
           x+=resultsRecovery[k][l];
           X[1].push_back(x);
           l++;
        }
        else {
            Y[1].push_back(coeff[k][i]);
            Y[1].push_back(coeff[k][i]);
            X[1].push_back(x);
            x+=resultsFailure[k][m];
            X[1].push_back(x);
            m++;
        }
    }
    X[1].push_back(x+resultsFailure[k][resultsFailure[k].size()-1]);
    Y[1].push_back(coeff[k][coeff[k].size()-1]);

    double sumt = 0;
    for(int i = 0; i < resultsFailure[k].size(); i++) {
        sumt += resultsFailure[k][i];
    }
    for(int i = 0; i < resultsRecovery[k].size(); i++) {
        sumt += resultsRecovery[k][i];
    }

}
//График без отказной работы
void GrafickWithoutRecovery(QVector<QVector<double>>& coeff,int k) {
    Y[2].push_back(coeff[k][0]*1.);
    for(int i = 1; i < coeff[k].size(); i++) {
        if(i%2 == 0) {
            Y[2].push_back(coeff[k][i]);
        }
    }
    double x = 0;
    double step = sim.TimeSimulation / Y[2].size();
    for(int i = 0; i < Y[2].size(); i++) {
        X[2].push_back(x);
        x += step;
    }
}

//расчет аналитического коэффиццинтеа
double TotalCoeffSystemAnalytical(QVector<QVector<double>> MatrixElements, QVector<QVector<double>> &CoeffsOfElement, QVector<QVector<double>> &CoeffsForEachElements, QVector<QVector<double>> &resultFailure, QVector<QVector<double>> &resultRecovery, QVector<QVector<double>>& BinaryElements) {
    int k = 0;
    for(int i = 0; i < MatrixElements.size(); i++) {
        CoeffsForEachElements.push_back(QVector<double>());
        for(int j = 0; j < MatrixElements[i].size(); j++) {
            if(MatrixElements[i][j] == 1) {
                CoeffsOfElement.push_back(QVector<double>());
                resultRecovery.push_back(QVector<double>());
                resultFailure.push_back(QVector<double>());

                BinaryElements.push_back(QVector<double>());
                CoeffsForEachElements[i].push_back(calculateCoeffOfElement(CoeffsOfElement[k],resultRecovery[k],resultFailure[k],BinaryElements[k]));
                k++;
            }
            else {
                CoeffsForEachElements[i].push_back(1);
            }
        }
    }
    QVector<double> temp;
    double pr = 1;
    //считаем коэффициент по столбцам
    for(int i = 0; i < CoeffsForEachElements.size(); i++) {
        for(int j = 0; j < CoeffsForEachElements[i].size(); j++){
            pr*=CoeffsForEachElements[i][j];
        }
        temp.push_back(1 - pr);
        pr = 1;
    }
    pr = 1;
    //считаем по строкам
    for(int i = 0; i < temp.size(); i++) {
        pr*= temp[i];
    }
    return 1 - pr;
}




void MainWindow::on_pushButton_2_clicked() {
    srand(time(NULL));
    QVector<QVector<double>> MatrixElements;
    //можно вывести график коэффициента для любого элемента i - элемент, j - коэффициенты
    QVector<QVector<double>>CoeffsOfElementGraph;
   QVector<QVector<double>>CoeffsOf;
   QVector<QVector<double>> CoeffsAnalytics;
    //матрица коэффициентов для каждого элемента
    QVector<QVector<double>> MatrixCoeffElements;
    QVector<QVector<double>> resultRecoveryOfElements;
    QVector<QVector<double>> resultFailureOfElements;
    QVector<QVector<int>> W;
    QVector<double> KrForWay;
    vector<pair<int,int>> WaysPair;

    vector<vector<int>>().swap(graph);
    vector<vector<int>>().swap(routs);
    vector<int>().swap(path);
    set<set<int>>().swap(ribs);
    ui->plainTextEdit->clear();

    QString inputText= ui->lineEdit_7->text();
    readFile(inputText, W);
    if(isfile) {
        ui->plainTextEdit->setPlainText(ui->plainTextEdit->toPlainText() + "Неверное имя файла!\n");
        return;
    }
    else if(W.size() == 0){
        ui->plainTextEdit->setPlainText(ui->plainTextEdit->toPlainText() + "Файл пуст, либо в нем некорретные данные!\n");
        return;
    }
    //создаем структуру из матрицы смежности
    for (int i = 0; i < sqrt(W.size()); i++) {
        MatrixElements.push_back(QVector<double>());
        for (int j = 0; j < sqrt(W[i].size()); j++) {
            MatrixElements[i].push_back(1);
        }
    }
    QString str = printStruct(MatrixElements);
    ui->plainTextEdit_3->setPlainText(ui->plainTextEdit_3->toPlainText() + str);
    if(ui->lineEdit_2->text().isEmpty() || ui->lineEdit_3->text().isEmpty()){
        ui->plainTextEdit->setPlainText(ui->plainTextEdit->toPlainText() + "Ошибка! Введите все параметры восстановления. \n");
        return;
    }
    else {
        sim.muRecovery = ui->lineEdit_2->text().toFloat();
        sim.sigmaRecovery = ui->lineEdit_3->text().toFloat();
    }
    if(ui->lineEdit_4->text().isEmpty() || ui->lineEdit_5->text().isEmpty()){
        ui->plainTextEdit->setPlainText(ui->plainTextEdit->toPlainText() + "Ошибка! Введите все параметры отказа. \n");
        return;
    }
    else {
        sim.muFailure = ui->lineEdit_4->text().toFloat();
        sim.sigmaFailure = ui->lineEdit_5->text().toFloat();
    }
    if(ui->lineEdit_6->text().isEmpty()){
        ui->plainTextEdit->setPlainText(ui->plainTextEdit->toPlainText() + "Ошибка! Введите количество итераций. \n");
        return;
    }
    else {
        sim.countIterations = ui->lineEdit_6->text().toFloat();
    }
    if(ui->lineEdit_9->text().isEmpty()){
        ui->plainTextEdit->setPlainText(ui->plainTextEdit->toPlainText() + "Ошибка! Введите время симуляции. \n");
        return;
    }
    else {
        sim.TimeSimulation = ui->lineEdit_9->text().toFloat();
    }

    QVector<QVector<double>>().swap(X);
    QVector<QVector<double>>().swap(Y);

    X.resize(8);
    Y.resize(8);
    QVector<QVector<double>> BinaryElements;
    QVector<int> ways;
    // находим все пути и выбираем кратчайшие
    spisSmej(W);
    searchWays(MatrixElements);
    int zn=0;
    if(MatrixElements.size()==2)
    {
        zn=0;
    }
    else{
        zn=MatrixElements.size()-1;
    }
    checkWays(ways, routs,zn);
    cout << "programm work!" << endl;
    double kr;
    //расчитываем аналитический коэффициент
    kr = TotalCoeffSystemAnalytical(MatrixElements,CoeffsOfElementGraph,MatrixCoeffElements,resultFailureOfElements,resultRecoveryOfElements,BinaryElements);
    ui->plainTextEdit_2->setPlainText("Аналитический коэффициент готовности " + QString::number(kr) + "\n");

//расчет иммитационных коэффов

    double x = TotalCoeffSystemImitationLower(BinaryElements, ways,KrForWay,MatrixElements);

    double y = TotalCoeffSystemImitationUpper(BinaryElements,ways,MatrixElements);

    cout << "programm work!" << endl;
   string out = "Кратшайшие пути \n";
    for(int i = 0; i < ways.size(); i++) {
        for (int j = 0;j<routs[ways[i]].size() - 1; j++) {
             out += "(" + to_string(routs[ways[i]][j]) + "," + to_string(routs[ways[i]][j + 1]) + ")";
        }
        out += " Kr way = " + to_string(KrForWay[i])+"\n";
    }
    QString ent = QString::fromStdString(out);
    ui->plainTextEdit->setPlainText(ui->plainTextEdit->toPlainText() + ent);

    ui->plainTextEdit_2->setPlainText(ui->plainTextEdit_2->toPlainText() + "Имитационный коэффициент готовности нижний  " + QString::number(x) + "\n");
    ui->plainTextEdit_2->setPlainText(ui->plainTextEdit_2->toPlainText() + "Имитационный коэффициент готовности верхний   " + QString::number(y) + "\n");
    int k = 1;
    cout << "Coeff for imit Lower=" << x << endl;
    cout << "Coeff for imit Upper=" << y << endl;
    cout << "Coeff for imit Sred=" << (x + y) / 2 << endl;

    CoeffsOf.resize(CoeffsOfElementGraph.size());

    for(int i = 0; i < Y[3].size(); i++) {
        Y[5].push_back((Y[3][i] + Y[4][i]) / 2);
        X[5].push_back(X[3][i]);
    }



    cout<<"I size"<<CoeffsOfElementGraph[k].size()<<endl;
    //CoeffsOf[k].resize(CoeffsOfElementGraph[k].size());

    CoeffAnalytic(CoeffsOf[k],CoeffsOfElementGraph[k].size());
    cout<<"А size"<<CoeffsOf[k].size()<<endl;
    //printVD("CoeffsA",CoeffsOf[k]);
   // GrafickCoeff(CoeffsOfElementGraph[k],sim.TimeSimulation);
  // double step = sim.TimeSimulation/(double) CoeffsOfElementGraph[k].size();
   /* double x0=0;
    Y[0].push_back(1);
    X[0].push_back(0);

    x0+=sim.muFailure;
     for(int i=1;i<CoeffsOfElementGraph[k].size()-1;i++)
     {

         if(CoeffsOfElementGraph[k][i]>CoeffsOfElementGraph[k][i+1])
         {
             x0+=sim.muFailure;
            Y[0].push_back(CoeffsOfElementGraph[k][i]);
             X[0].push_back(x0);
         }
         else {
             x0+=sim.muRecovery;
            Y[0].push_back(CoeffsOfElementGraph[k][i]);
             X[0].push_back(x0);

         }

     }

*/
    double b=0;
    double s=0;
    for (int i = CoeffsOfElementGraph[k].size()-1; i > CoeffsOfElementGraph[k].size()/2; i--) {
        b += CoeffsOfElementGraph[k][i];
        s++;
    }
    b = b / s;
    b -= 0.0001;



    Y[0].push_back(1);
    double x0=0;
    X[0].push_back(x0);
    x0+=sim.muFailure;
 /*   for (int i = availabilityFactor.Count - 2; i > 0; i--)
               {
                   if (Math.Abs(b - availabilityFactor[i]) >= 0.04)
                   {
                       goto Leave;
                   }
                   else
                   {
                       if (Math.Abs(b - availabilityFactor[i]) >= 0.003)
                       {
                           break;
                       }
                       else
                       {
                           goto Leave;
                       }

                   }
               }*/
   /* for(int i=1;i<CoeffsOfElementGraph[k].size()-1;i++)
    {
    if(CoeffsOfElementGraph[k][i]>CoeffsOfElementGraph[k][i+1])
    {
        if(abs(CoeffsOfElementGraph[k][i]-CoeffsOfElementGraph[k][i+1]) >= 0.04)
        {
        //x0+=sim.muRecovery;
            Y[0].push_back(CoeffsOfElementGraph[k][i]);
             X[0].push_back(x0);
             x0+=sim.muFailure;
        }
        else{
            if(abs(b-CoeffsOfElementGraph[k][i] >= 0.003))
            {
                Y[0].push_back(CoeffsOfElementGraph[k][CoeffsOfElementGraph[k].size()-1]);
                 X[0].push_back(x0);
                 x0+=sim.muFailure;
                }
            }
            /*else{
            Y[0].push_back(CoeffsOfElementGraph[k][CoeffsOfElementGraph[k].size()-1]);
             X[0].push_back(x0);
             x0+=sim.muFailure;
            }
        }

    else {

            if(abs(b-CoeffsOfElementGraph[k][i] >= 0.003))
            {
                Y[0].push_back(CoeffsOfElementGraph[k][CoeffsOfElementGraph[k].size()-1]);
                 X[0].push_back(x0);
                 x0+=sim.muRecovery;
            }
            /*else{
            Y[0].push_back(CoeffsOfElementGraph[k][CoeffsOfElementGraph[k].size()-1]);
             X[0].push_back(x0);
             x0+=sim.muRecovery;
            }
        }

    }*/
    for(int i=1;i<CoeffsOfElementGraph[k].size()/2;i++)
        {
        if(CoeffsOfElementGraph[k][i]>CoeffsOfElementGraph[k][i+1])
        {
           // if(abs(CoeffsOfElementGraph[k][i]-CoeffsOfElementGraph[k][i+1])>=0.04){
            //x0+=sim.muRecovery;
            Y[0].push_back(CoeffsOfElementGraph[k][i]);
            X[0].push_back(x0);
            x0+=sim.muFailure;
         /*   }
            else{
                Y[0].push_back(CoeffsOfElementGraph[k][CoeffsOfElementGraph[k].size()-1]);
                 X[0].push_back(x0);
                 x0+=sim.muFailure;
            }*/
        }
        else {
            //if(abs(b-CoeffsOfElementGraph[k][i] >= 0.04))
           // {
                Y[0].push_back(CoeffsOfElementGraph[k][i]);
                 X[0].push_back(x0);
                 x0+=sim.muRecovery;
           // }
           /* else{
            Y[0].push_back(CoeffsOfElementGraph[k][CoeffsOfElementGraph[k].size()-1]);
             X[0].push_back(x0);
             x0+=sim.muRecovery;
            }
        }*/
        }
    }
    for(int i=CoeffsOfElementGraph[k].size()/2;i<CoeffsOfElementGraph[k].size()-1;i++)
    {
        if(CoeffsOfElementGraph[k][i]>CoeffsOfElementGraph[k][i+1])
        {
           // if(abs(CoeffsOfElementGraph[k][i]-CoeffsOfElementGraph[k][i+1])>=0.04){
            //x0+=sim.muRecovery;
            Y[0].push_back(b);
            X[0].push_back(x0);
            x0+=sim.muFailure;
         /*   }
            else{
                Y[0].push_back(CoeffsOfElementGraph[k][CoeffsOfElementGraph[k].size()-1]);
                 X[0].push_back(x0);
                 x0+=sim.muFailure;
            }*/
        }
        else {

                Y[0].push_back(b);
                 X[0].push_back(x0);
                 x0+=sim.muRecovery;

        }
    }
   printVD("Coeff A betwin erase",CoeffsOf[k]);
   CoeffsOf[k].erase(CoeffsOf[k].begin()+1);
   printVD("Coeff A after erase",CoeffsOf[k]);
     double step2 = sim.TimeSimulation/(double) CoeffsOf[k].size();
     double x6=0;
     Y[6].push_back(1);
     X[6].push_back(x6);
     x6+=sim.muFailure;

     /* for(int i=0;i<CoeffsOf[k].size();i++)
      {
          Y[6].push_back(CoeffsOf[k][i]);
          X[6].push_back(x6);
          x6+=step2;
      }*/
      for(int i=1;i<CoeffsOf[k].size()-1;i++)
      {

          if(CoeffsOf[k][i]>CoeffsOf[k][i+1])
          {

             Y[6].push_back(CoeffsOf[k][i]);
              X[6].push_back(x6);
              x6+=sim.muFailure;
          }
          else {

             Y[6].push_back(CoeffsOf[k][i]);
              X[6].push_back(x6);
              x6+=sim.muRecovery;

          }

      }
      Y[6].push_back(CoeffsOf[k][CoeffsOf[k].size()-1]);
      Y[6].push_back(CoeffsOf[k][CoeffsOf[k].size()-1]);
      X[6].push_back(x6+sim.muFailure);
      x6=x6+sim.muFailure;
      X[6].push_back(x6+sim.muFailure);
      x6=x6+sim.muFailure;
     Y[0].push_back(b);
     Y[0].push_back(b);
     X[0].push_back(x0+sim.muFailure);
      x0+=sim.muFailure;
     X[0].push_back(x0+sim.muFailure);
      Y[6].push_back(CoeffsOf[k][CoeffsOf[k].size()-1]);
      Y[6].push_back(CoeffsOf[k][CoeffsOf[k].size()-1]);
      X[6].push_back(x6+sim.muFailure);
      x6=x6+sim.muFailure;
      X[6].push_back(x6+sim.muFailure);

    /*  QVector<double> x1;
      QVector<double> y2;

      interpol(x1,y2);

      X[6].swap(x1);
      Y[6].swap(y2);

      QVector<double> x12;
      QVector<double> y22;

     interpol1(x12,y22);

     X[0].swap(x12);
     Y[0].swap(y22);
*/


   // GrafickCoeffA(CoeffsOf[k],sim.TimeSimulation);
    printVD("CoeffsA",CoeffsOf[k]);
    // GrafickCoeff(CoeffsOfElementGraph[k],sim.TimeSimulation);
     printVD("CoeffsI",CoeffsOfElementGraph[k]);
     cout<<"Y0size"<<Y[0].size()<<endl;
     cout<<"Y6size"<<Y[6].size()<<endl;

      printVD("CoeffsI",Y[0]);
      cout<<"B = "<<b<<endl;
    GrafickTime(resultRecoveryOfElements,resultFailureOfElements,CoeffsOfElementGraph,k);
    GrafickWithoutRecovery(CoeffsOfElementGraph,k);
    cout<<"wordk end!"<<endl;
}

void MainWindow::on_pushButton_4_clicked() {
    //ui->plainTextEdit->clear();
    if(X.size() == 0) {
        ui->plainTextEdit->setPlainText("Нажмите сначала на кнопку рассчитать коэффициент!");
        return;
    }
    double a = X[1][0]; //Начало интервала, где рисуем график по оси Ox
    double b = X[1].back();
    ui->widget->clearGraphs();//Если нужно, но очищаем все графики
    //Добавляем один график в widget
    ui->widget->addGraph();
    //Говорим, что отрисовать нужно график по нашим двум массивам x и y
    ui->widget->graph(0)->setData(X[1], Y[1]);
    //Подписываем оси Ox и Oy
    ui->widget->xAxis->setLabel("Время в часах");
    ui->widget->yAxis->setLabel("Коэфф. готовности ");
    //Установим область, которая будет показываться на графике
    ui->widget->xAxis->setRange(a, b);//Для оси Ox
    //Для показа границ по оси Oy сложнее, так как надо по правильному
    //вычислить минимальное и максимальное значение в векторах
    ui->widget->yAxis->setRange(Y[1][1]-0.0001, 1);//Для оси Oy
    //И перерисуем график на нашем widget
    ui->widget->replot();
}


void MainWindow::on_pushButton_3_clicked() {
    //ui->plainTextEdit->clear();
    if(X.size() == 0) {
        ui->plainTextEdit->setPlainText("Нажмите сначала на кнопку рассчитать коэффициент!");
        return;
    }
    double a = X[0][0]; //Начало интервала, где рисуем график по оси Ox
    double b = X[0].back();
    ui->widget->clearGraphs();//Если нужно, но очищаем все графики
    //Добавляем один график в widget
    ui->widget->addGraph();
    ui->widget->addGraph();
    ui->widget->graph(1)->setPen(QPen(Qt::red));


    //ui->widget->graph(0)->setAntialiased(QCP::AntialiasedElement(aeAll));
    //Говорим, что отрисовать нужно график по нашим двум массивам x и y
    ui->widget->graph(0)->setData(X[0], Y[0]);
    ui->widget->graph(1)->setData(X[6], Y[6]);
    double min=1000;
    for(int i=0;i<Y[6].size();i++)
    {
        if(Y[6][i]<min)
        {
            min=Y[6][i];
        }
    }
    for(int i=0;i<Y[0].size();i++)
    {
        if(Y[0][i]<min)
        {
            min=Y[0][i];
        }
    }

    //Подписываем оси Ox и Oy
    ui->widget->xAxis->setLabel("Время в часах");
    ui->widget->yAxis->setLabel("Коэффициент готовности ");
    //Установим область, которая будет показываться на графике
    ui->widget->xAxis->setRange(a, b);//Для оси Ox
    //Для показа границ по оси Oy сложнее, так как надо по правильному
    //вычислить минимальное и максимальное значение в векторах

    ui->widget->yAxis->setRange(min-0.01, 1.1);;//Для оси Oy
    //И перерисуем график на нашем widget
    ui->widget->replot();
}

void MainWindow::on_pushButton_6_clicked() {
    //ui->plainTextEdit->clear();
    if(X.size() == 0) {
        ui->plainTextEdit->setPlainText("Нажмите сначала на кнопку рассчитать коэффициент!");
        return;
    }
    double a = X[3][0]; //Начало интервала, где рисуем график по оси Ox
    double b = X[3].back();
    ui->widget->clearGraphs();//Если нужно, но очищаем все графики
    //Добавляем один график в widget
    ui->widget->addGraph();
    //Говорим, что отрисовать нужно график по нашим двум массивам x и y
    ui->widget->graph(0)->setData(X[3], Y[3]);
    ui->widget->addGraph();
    //Говорим, что отрисовать нужно график по нашим двум массивам x и y
    ui->widget->graph(1)->setData(X[4], Y[4]);
    ui->widget->addGraph();
    //Говорим, что отрисовать нужно график по нашим двум массивам x и y
    ui->widget->graph(2)->setData(X[5], Y[5]);
    ui->widget->graph(0)->setPen(QPen(Qt::red));
    ui->widget->graph(1)->setPen(QPen(Qt::blue));
    ui->widget->graph(2)->setPen(QPen(Qt::black));
    ui->plainTextEdit_2->setPlainText(ui->plainTextEdit_2->toPlainText() + "Нижняя граница - Синий "+ "\n");
    ui->plainTextEdit_2->setPlainText(ui->plainTextEdit_2->toPlainText() + "Верхняя граница - Красный"+ "\n");
    ui->plainTextEdit_2->setPlainText(ui->plainTextEdit_2->toPlainText() + "Среднее - Черный "+ "\n");
    //Подписываем оси Ox и Oy
    ui->widget->xAxis->setLabel("Время в минутах");
    ui->widget->yAxis->setLabel("Вероятность безотказной работы ");
    //Установим область, которая будет показываться на графике
    ui->widget->xAxis->setRange(a, b);//Для оси Ox
    //Для показа границ по оси Oy сложнее, так как надо по правильному
    //вычислить минимальное и максимальное значение в векторах
    ui->widget->yAxis->setRange(Y[3][0], 1);//Для оси Oy
    //И перерисуем график на нашем widget
    ui->widget->replot();
}








void MainWindow::on_pushButton_5_clicked()
{
   // ui->plainTextEdit->clear();
    if(X.size() == 0) {
        ui->plainTextEdit->setPlainText("Нажмите сначала на кнопку рассчитать коэффициент!");
        return;
    }
    double a = X[6][0]; //Начало интервала, где рисуем график по оси Ox
    double b = X[6].back();
    ui->widget->clearGraphs();//Если нужно, но очищаем все графики
    //Добавляем один график в widget
    ui->widget->addGraph();

    ui->widget->graph(0)->setAntialiased(true);
    //Говорим, что отрисовать нужно график по нашим двум массивам x и y
    ui->widget->graph(0)->setData(X[6], Y[6]);
    ui->widget->graph(0)->setAntialiased(true);
    cout<<"Y0 size = "<<Y[6].size()<<endl;
    printVD("Y0",Y[6]);
    printVD("X0",X[6]);
    //Подписываем оси Ox и Oy
    ui->widget->xAxis->setLabel("Время в часах");
    ui->widget->yAxis->setLabel("Вероятность безотказной работы ");
    //Установим область, которая будет показываться на графике
    ui->widget->xAxis->setRange(a, b);//Для оси Ox
    //Для показа границ по оси Oy сложнее, так как надо по правильному
    //вычислить минимальное и максимальное значение в векторах
    ui->widget->yAxis->setRange(Y[6][Y[6].size()-1]-0.001, 1);//Для оси Oy
    ui->widget->graph(0)->setAntialiased(true);
    //И перерисуем график на нашем widget
    ui->widget->replot();
}
