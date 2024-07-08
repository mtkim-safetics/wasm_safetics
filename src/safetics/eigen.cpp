#include <emscripten/bind.h>
#include <emscripten/val.h>

#include <sstream>

#define DECOMPOSITIONS

#include "../external/Eigen/Core"
#include "../external/Eigen/LU"

#ifdef DECOMPOSITIONS
#include "../external/Eigen/Cholesky"
#include "../external/Eigen/Eigenvalues"
#include "../external/Eigen/QR"
#include "../external/Eigen/SVD"
#endif// DECOMPOSITIONS

using namespace emscripten;

typedef Eigen::MatrixXf MatrixXf;
typedef Eigen::Diagonal<Eigen::MatrixXf, Eigen::DynamicIndex> DiagonalXf;
typedef Eigen::Block<MatrixXf, Eigen::Dynamic, Eigen::Dynamic, false> ColMajorBlockXf;
typedef Eigen::Block<MatrixXf, Eigen::Dynamic, Eigen::Dynamic, true> RowMajorBlockXf;
typedef Eigen::Block<MatrixXf, Eigen::Dynamic, 1, true> ColXf;
typedef Eigen::Block<MatrixXf, 1, Eigen::Dynamic, false> RowXf;

#ifdef DECOMPOSITIONS
typedef Eigen::FullPivLU<Eigen::MatrixXf> FullPivLUXf;
typedef Eigen::PartialPivLU<Eigen::MatrixXf> PartialPivLUXf;
typedef Eigen::HouseholderQR<Eigen::MatrixXf> HouseholderQRXf;
typedef Eigen::ColPivHouseholderQR<Eigen::MatrixXf> ColPivHouseholderQRXf;
typedef Eigen::FullPivHouseholderQR<Eigen::MatrixXf> FullPivHouseholderQRXf;
typedef Eigen::BDCSVD<Eigen::MatrixXf> BDCSVDXf;
typedef Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> SelfAdjointEigenSolverXf;
#endif// DECOMPOSITIONS

struct Shape {
    int rows, cols;
};

EMSCRIPTEN_BINDINGS(EigenJS) {
    value_array<Shape>("Shape")
            .element(&Shape::rows)
            .element(&Shape::cols);

    constant("ComputeFullU", (int) Eigen::DecompositionOptions::ComputeFullU);
    constant("ComputeThinU", (int) Eigen::DecompositionOptions::ComputeThinU);
    constant("ComputeFullV", (int) Eigen::DecompositionOptions::ComputeFullV);
    constant("ComputeThinV", (int) Eigen::DecompositionOptions::ComputeThinV);
    constant("ComputeEigenvectors", (int) Eigen::DecompositionOptions::ComputeEigenvectors);
    constant("EigenvaluesOnly", (int) Eigen::DecompositionOptions::EigenvaluesOnly);

    class_<MatrixXf>("Matrixf")
            .constructor<int, int>()
            .class_function("fromArray", optional_override([](const emscripten::val& v, Shape shape) {
                                // This function uses the technique from this issue but without first passing through a std::vector.
                                //    https://github.com/emscripten-core/emscripten/pull/5655#issuecomment-520378179
                                const auto l = v["length"].as<unsigned>();
                                assert(l == shape.rows * shape.cols);
                                MatrixXf* m = new MatrixXf(shape.rows, shape.cols);
                                emscripten::val memoryView{emscripten::typed_memory_view(l, m->data())};
                                memoryView.call<void>("set", v);
                                return m;
                            }),
                            allow_raw_pointers())
            .class_function("columnVector", optional_override([](const emscripten::val& v) {
                                // This function uses the technique from this issue but without first passing through a std::vector.
                                //    https://github.com/emscripten-core/emscripten/pull/5655#issuecomment-520378179
                                const auto l = v["length"].as<unsigned>();
                                MatrixXf* m = new MatrixXf(l, 1);
                                emscripten::val memoryView{emscripten::typed_memory_view(l, m->data())};
                                memoryView.call<void>("set", v);
                                return m;
                            }),
                            allow_raw_pointers())
            .class_function("rowVector", optional_override([](const emscripten::val& v) {
                                // This function uses the technique from this issue but without first passing through a std::vector.
                                //    https://github.com/emscripten-core/emscripten/pull/5655#issuecomment-520378179
                                const auto l = v["length"].as<unsigned>();
                                MatrixXf* m = new MatrixXf(1, l);
                                emscripten::val memoryView{emscripten::typed_memory_view(l, m->data())};
                                memoryView.call<void>("set", v);
                                return m;
                            }),
                            allow_raw_pointers())
            .class_function("from", optional_override([](emscripten::val valueFn, Shape dims) {
                                int m = dims.rows;
                                int n = dims.cols;
                                MatrixXf* a = new MatrixXf(m, n);
                                for (int j = 0; j < n; j++) {
                                    for (int i = 0; i < m; i++) {
                                        (*a)(i, j) = valueFn(i, j).as<float>();
                                    }
                                }
                                return a;
                            }),
                            allow_raw_pointers())
            .class_function("fromDiagonal", optional_override([](const MatrixXf& vector) {
                                if (vector.cols() == 1) {
                                    return new MatrixXf(vector.asDiagonal());
                                } else {
                                    return new MatrixXf(vector.transpose().asDiagonal());
                                }
                            }),
                            allow_raw_pointers())
            .class_function("fromBlock", optional_override([](const ColMajorBlockXf& block) {
                                return new MatrixXf(block);
                            }),
                            allow_raw_pointers())
            .class_function("fromRowMajorBlock", optional_override([](const RowMajorBlockXf& block) {
                                return new MatrixXf(block);
                            }),
                            allow_raw_pointers())
            .class_function("fromCol", optional_override([](const ColXf& col) {
                                return new MatrixXf(col);
                            }),
                            allow_raw_pointers())
            .class_function("fromRow", optional_override([](const RowXf& row) {
                                return new MatrixXf(row);
                            }),
                            allow_raw_pointers())
            .class_function("random", optional_override([](Shape dims) {
                                MatrixXf A(MatrixXf::Random(dims.rows, dims.cols));
                                return A;
                            }))
            .class_function("zeros", optional_override([](Shape dims) {
                                MatrixXf A(MatrixXf::Zero(dims.rows, dims.cols));
                                return A;
                            }))
            .class_function("ones", optional_override([](Shape dims) {
                                MatrixXf A(MatrixXf::Ones(dims.rows, dims.cols));
                                return A;
                            }))
            .class_function("identity", optional_override([](Shape dims) {
                                MatrixXf A(MatrixXf::Identity(dims.rows, dims.cols));
                                return A;
                            }))

            .property("rows", optional_override([](const MatrixXf& self) { return self.rows(); }))
            .property("cols", optional_override([](const MatrixXf& self) { return self.cols(); }))
            .property("shape", optional_override([](const MatrixXf& self) {
                          return Shape{self.rows(), self.cols()};
                      }))
            .property("buffer", optional_override([](const MatrixXf& self) {
                          return val(typed_memory_view(self.size(), (float*) self.data()));
                      }))
            .property("offset", optional_override([](const MatrixXf& self) {
                          return 0;
                      }))
            .property("outerStride", optional_override([](const MatrixXf& self) {
                          return self.outerStride();
                      }))
            .property("innerStride", optional_override([](const MatrixXf& self) {
                          return self.innerStride();
                      }))
            .property("stride", optional_override([](const MatrixXf& self) {
                          return Shape{self.innerStride(), self.outerStride()};
                      }))
            .property("dtype", optional_override([](const MatrixXf& self) {
                          return std::string("float32");
                      }))

            .function("block", optional_override([](MatrixXf& self, int startRow, int startCol, Shape dims) {
                          return new Eigen::Block<MatrixXf, Eigen::Dynamic, Eigen::Dynamic>(self, startRow, startCol, dims.rows, dims.cols);
                      }),
                      allow_raw_pointers())
            .function("col", optional_override([](MatrixXf& self, int i) {
                          return self.col(i);
                      }),
                      allow_raw_pointers())
            .function("row", optional_override([](MatrixXf& self, int i) {
                          return self.row(i);
                      }),
                      allow_raw_pointers())

            .function("topLeftCorner", optional_override([](MatrixXf& self, Shape dims) {
                          return new ColMajorBlockXf(self.topLeftCorner(dims.rows, dims.cols));
                      }),
                      allow_raw_pointers())
            .function("topRightCorner", optional_override([](MatrixXf& self, Shape dims) {
                          return new ColMajorBlockXf(self.topRightCorner(dims.rows, dims.cols));
                      }),
                      allow_raw_pointers())
            .function("bottomLeftCorner", optional_override([](MatrixXf& self, Shape dims) {
                          return new ColMajorBlockXf(self.bottomLeftCorner(dims.rows, dims.cols));
                      }),
                      allow_raw_pointers())
            .function("bottomRightCorner", optional_override([](MatrixXf& self, Shape dims) {
                          return new ColMajorBlockXf(self.bottomRightCorner(dims.rows, dims.cols));
                      }),
                      allow_raw_pointers())

            .function("topRows", optional_override([](MatrixXf& self, int n) {
                          return new ColMajorBlockXf(self.topRows(n));
                      }),
                      allow_raw_pointers())
            .function("bottomRows", optional_override([](MatrixXf& self, int n) {
                          return new ColMajorBlockXf(self.bottomRows(n));
                      }),
                      allow_raw_pointers())
            .function("leftCols", optional_override([](MatrixXf& self, int n) {
                          return new RowMajorBlockXf(self.leftCols(n));
                      }),
                      allow_raw_pointers())
            .function("rightCols", optional_override([](MatrixXf& self, int n) {
                          return new RowMajorBlockXf(self.rightCols(n));
                      }),
                      allow_raw_pointers())

            .function("setConstant", optional_override([](MatrixXf& self, float value) {
                          self.setConstant(value);
                          return self;
                      }))
            .function("setZero", optional_override([](MatrixXf& self) {
                          self.setZero();
                          return self;
                      }))
            .function("setOnes", optional_override([](MatrixXf& self) {
                          self.setOnes();
                          return self;
                      }))
            .function("setIdentity", optional_override([](MatrixXf& self) {
                          self.setIdentity();
                          return self;
                      }))
            .function("setRandom", optional_override([](MatrixXf& self) {
                          self.setRandom();
                          return self;
                      }))

            .function("mutate", optional_override([](MatrixXf& self, emscripten::val mutateFn) {
                          int m = self.rows();
                          int n = self.cols();
                          for (int j = 0; j < n; j++) {
                              for (int i = 0; i < m; i++) {
                                  self(i, j) = mutateFn(i, j, self(i, j)).as<float>();
                              }
                          }
                          return &self;
                      }),
                      allow_raw_pointers())
            .function("map", optional_override([](const MatrixXf& self, emscripten::val valueFn) {
                          int m = self.rows();
                          int n = self.cols();
                          MatrixXf* a = new MatrixXf(m, n);
                          for (int j = 0; j < n; j++) {
                              for (int i = 0; i < m; i++) {
                                  (*a)(i, j) = valueFn(i, j, self(i, j)).as<float>();
                              }
                          }
                          return a;
                      }),
                      allow_raw_pointers())
            .function("transposeInPlace", optional_override([](MatrixXf& self) {
                          self.transposeInPlace();
                          return &self;
                      }),
                      allow_raw_pointers())
            .function("transpose", optional_override([](MatrixXf& self) {
                          MatrixXf transposed(self.transpose());
                          return transposed;
                      }))
            .function("diagonal", optional_override([](MatrixXf& self, int offset) {
                          return new DiagonalXf(self.diagonal(offset));
                      }),
                      allow_raw_pointers())
            .function("get", optional_override([](const MatrixXf& self, int i, int j) {
                          return self(i, j);
                      }))
            .function("set", optional_override([](MatrixXf& self, int i, int j, float value) {
                          self(i, j) = value;
                          return self;
                      }))
            .function("determinant", optional_override([](const MatrixXf& self) {
                          return self.determinant();
                      }))
            .function("inverse", optional_override([](const MatrixXf& self) {
                          MatrixXf* m = new MatrixXf(self.inverse());
                          return m;
                      }),
                      allow_raw_pointers())
            .function("squaredNorm", optional_override([](const MatrixXf& self) {
                          return self.squaredNorm();
                      }))
            .function("norm", optional_override([](const MatrixXf& self) {
                          return self.norm();
                      }))
            .function("l1Norm", optional_override([](const MatrixXf& self) {
                          return self.lpNorm<1>();
                      }))
            .function("lInfNorm", optional_override([](const MatrixXf& self) {
                          return self.lpNorm<Eigen::Infinity>();
                      }))

#ifdef DECOMPOSITIONS
            .function("fullPivLU", optional_override([](MatrixXf& self) {
                          return new FullPivLUXf(self);
                      }),
                      allow_raw_pointers())
            .function("ldlt", optional_override([](MatrixXf& self) {
                          return new Eigen::LDLT<MatrixXf>(self);
                      }),
                      allow_raw_pointers())
            .function("llt", optional_override([](MatrixXf& self) {
                          return new Eigen::LLT<MatrixXf>(self);
                      }),
                      allow_raw_pointers())
            .function("partialPivLU", optional_override([](MatrixXf& self) {
                          return new PartialPivLUXf(self);
                      }),
                      allow_raw_pointers())
            .function("householderQR", optional_override([](MatrixXf& self) {
                          return new HouseholderQRXf(self);
                      }),
                      allow_raw_pointers())
            .function("colPivHouseholderQR", optional_override([](MatrixXf& self) {
                          return new ColPivHouseholderQRXf(self);
                      }),
                      allow_raw_pointers())
            .function("fullPivHouseholderQR", optional_override([](MatrixXf& self) {
                          return new FullPivHouseholderQRXf(self);
                      }),
                      allow_raw_pointers())
            .function("bdcSVD", optional_override([](MatrixXf& self, unsigned int opts) {
                          return new BDCSVDXf(self, opts);
                      }),
                      allow_raw_pointers())
            .function("selfAdjointEigenSolver", optional_override([](MatrixXf& self, unsigned int opts) {
                          return new SelfAdjointEigenSolverXf(self, opts);
                      }),
                      allow_raw_pointers())

            // Shortcuts for preferred versions:
            .function("svd", optional_override([](MatrixXf& self, unsigned int opts) {
                          return new BDCSVDXf(self, opts);
                      }),
                      allow_raw_pointers())
            .function("lu", optional_override([](MatrixXf& self) {
                          return new PartialPivLUXf(self);
                      }),
                      allow_raw_pointers())
            .function("qr", optional_override([](MatrixXf& self) {
                          return new ColPivHouseholderQRXf(self);
                      }),
                      allow_raw_pointers())
#endif// DECOMPOSITIONS

            .function("mul", optional_override([](const MatrixXf& self, const MatrixXf rhs) {
                          return new MatrixXf(self * rhs);
                      }),
                      allow_raw_pointers())
            .function("muleq", optional_override([](MatrixXf& self, const MatrixXf rhs) {
                          self *= rhs;
                          return &self;
                      }),
                      allow_raw_pointers())
            .function("muls", optional_override([](const MatrixXf& self, float rhs) {
                          return new MatrixXf(self * rhs);
                      }),
                      allow_raw_pointers())
            .function("mulseq", optional_override([](MatrixXf& self, float rhs) {
                          self *= rhs;
                          return &self;
                      }),
                      allow_raw_pointers())

            .function("add", optional_override([](const MatrixXf& self, const MatrixXf rhs) {
                          return new MatrixXf(self + rhs);
                      }),
                      allow_raw_pointers())
            .function("addeq", optional_override([](MatrixXf& self, const MatrixXf rhs) {
                          self += rhs;
                          return &self;
                      }),
                      allow_raw_pointers())
            .function("adds", optional_override([](const MatrixXf& self, const float rhs) {
                          return new MatrixXf(self.array() + rhs);
                      }),
                      allow_raw_pointers())
            .function("addseq", optional_override([](MatrixXf& self, const float rhs) {
                          self.array() += rhs;
                          return &self;
                      }),
                      allow_raw_pointers())

            .function("sub", optional_override([](const MatrixXf& self, const MatrixXf rhs) {
                          return new MatrixXf(self - rhs);
                      }),
                      allow_raw_pointers())
            .function("subeq", optional_override([](MatrixXf& self, const MatrixXf rhs) {
                          self -= rhs;
                          return &self;
                      }),
                      allow_raw_pointers())
            .function("subs", optional_override([](const MatrixXf& self, float rhs) {
                          return new MatrixXf(self.array() - rhs);
                      }),
                      allow_raw_pointers())
            .function("subseq", optional_override([](MatrixXf& self, float rhs) {
                          self.array() -= rhs;
                          return &self;
                      }),
                      allow_raw_pointers())

            .function("toString", optional_override([](const MatrixXf& self) {
                          std::stringstream buf;
                          buf << self << std::endl;
                          return buf.str();
                      }));

    class_<DiagonalXf>("Diagonalf")
            .function("setConstant", optional_override([](DiagonalXf& self, float value) {
                          self.setConstant(value);
                          return self;
                      }))
            .function("setRandom", optional_override([](DiagonalXf& self) {
                          self.setRandom();
                          return self;
                      }))
            .function("setZero", optional_override([](DiagonalXf& self) {
                          self.setZero();
                          return self;
                      }))
            .function("setOnes", optional_override([](DiagonalXf& self) {
                          self.setOnes();
                          return self;
                      }))
            .function("set", optional_override([](DiagonalXf& self, int i, float value) {
                          self(i) = value;
                          return self;
                      }))
            .function("get", optional_override([](DiagonalXf& self, int i) {
                          return self(i);
                      }))
            .function("assign", optional_override([](DiagonalXf& self, const MatrixXf& src) {
                          self = src;
                          return self;
                      }));

    class_<ColMajorBlockXf>("ColMajorBlockf")
            .property("rows", optional_override([](const ColMajorBlockXf& self) { return self.rows(); }))
            .property("cols", optional_override([](const ColMajorBlockXf& self) { return self.cols(); }))
            .property("shape", optional_override([](const ColMajorBlockXf& self) {
                          return Shape{self.rows(), self.cols()};
                      }))
            .function("setConstant", optional_override([](ColMajorBlockXf& self, float value) {
                          self.setConstant(value);
                          return self;
                      }))
            .function("setRandom", optional_override([](ColMajorBlockXf& self) {
                          self.setRandom();
                          return self;
                      }))
            .function("setZero", optional_override([](ColMajorBlockXf& self) {
                          self.setZero();
                          return self;
                      }))
            .function("setOnes", optional_override([](ColMajorBlockXf& self) {
                          self.setOnes();
                          return self;
                      }))
            .function("setIdentity", optional_override([](ColMajorBlockXf& self) {
                          self.setIdentity();
                          return self;
                      }))
            .function("toString", optional_override([](const ColMajorBlockXf& self) {
                          std::stringstream buf;
                          buf << self << std::endl;
                          return buf.str();
                      }))
            .function("set", optional_override([](ColMajorBlockXf& self, int i, int j, float value) {
                          self(i, j) = value;
                          return self;
                      }))
            .function("get", optional_override([](ColMajorBlockXf& self, int i, int j) {
                          return self(i, j);
                      }))
            .function("assign", optional_override([](ColMajorBlockXf& self, const MatrixXf& src) {
                          self = src;
                          return self;
                      }));
    ;

    class_<RowMajorBlockXf>("RowMajorBlockf")
            .property("rows", optional_override([](const RowMajorBlockXf& self) { return self.rows(); }))
            .property("cols", optional_override([](const RowMajorBlockXf& self) { return self.cols(); }))
            .property("shape", optional_override([](const RowMajorBlockXf& self) {
                          return Shape{self.rows(), self.cols()};
                      }))
            .function("setConstant", optional_override([](RowMajorBlockXf& self, float value) {
                          self.setConstant(value);
                          return self;
                      }))
            .function("setRandom", optional_override([](RowMajorBlockXf& self) {
                          self.setRandom();
                          return self;
                      }))
            .function("setZero", optional_override([](RowMajorBlockXf& self) {
                          self.setZero();
                          return self;
                      }))
            .function("setOnes", optional_override([](RowMajorBlockXf& self) {
                          self.setOnes();
                          return self;
                      }))
            .function("setIdentity", optional_override([](RowMajorBlockXf& self) {
                          self.setIdentity();
                          return self;
                      }))
            .function("toString", optional_override([](const RowMajorBlockXf& self) {
                          std::stringstream buf;
                          buf << self << std::endl;
                          return buf.str();
                      }))
            .function("set", optional_override([](RowMajorBlockXf& self, int i, int j, float value) {
                          self(i, j) = value;
                          return self;
                      }))
            .function("set", optional_override([](RowMajorBlockXf& self, int i, int j) {
                          return self(i, j);
                      }))
            .function("assign", optional_override([](RowMajorBlockXf& self, MatrixXf& src) {
                          self = src;
                          return self;
                      }));
    ;

    class_<RowXf>("Rowf")
            .function("setConstant", optional_override([](RowXf& self, float value) {
                          self.setConstant(value);
                          return self;
                      }))
            .function("setRandom", optional_override([](RowXf& self) {
                          self.setRandom();
                          return self;
                      }))
            .function("setZero", optional_override([](RowXf& self) {
                          self.setZero();
                          return self;
                      }))
            .function("setOnes", optional_override([](RowXf& self) {
                          self.setOnes();
                          return self;
                      }))
            .function("set", optional_override([](RowXf& self, int i, float value) {
                          self(i) = value;
                          return self;
                      }))
            .function("get", optional_override([](RowXf& self, int i) {
                          return self(i);
                      }))
            .function("assign", optional_override([](RowXf& self, const MatrixXf& src) {
                          self = src;
                          return self;
                      }));

    class_<ColXf>("Colf")
            .function("setConstant", optional_override([](ColXf& self, float value) {
                          self.setConstant(value);
                          return self;
                      }))
            .function("setRandom", optional_override([](ColXf& self) {
                          self.setRandom();
                          return self;
                      }))
            .function("setZero", optional_override([](ColXf& self) {
                          self.setZero();
                          return self;
                      }))
            .function("setOnes", optional_override([](ColXf& self) {
                          self.setOnes();
                          return self;
                      }))
            .function("set", optional_override([](ColXf& self, int i, float value) {
                          self(i) = value;
                          return self;
                      }))
            .function("get", optional_override([](ColXf& self, int i) {
                          return self(i);
                      }))
            .function("assign", optional_override([](ColXf& self, const MatrixXf& src) {
                          self = src;
                          return self;
                      }));

#ifdef DECOMPOSITIONS
    class_<Eigen::LDLT<MatrixXf>>("LDLTf")
            .function("matrixL", optional_override([](const Eigen::LDLT<MatrixXf>& self) {
                          return new MatrixXf(self.matrixL());
                      }),
                      allow_raw_pointers())
            .function("matrixU", optional_override([](const Eigen::LDLT<MatrixXf>& self) {
                          return new MatrixXf(self.matrixU());
                      }),
                      allow_raw_pointers())
            .function("transpositionsP", optional_override([](const Eigen::LDLT<MatrixXf>& self) {
                          MatrixXf* m = new MatrixXf(self.rows(), self.rows());
                          m->setIdentity();
                          *m = self.transpositionsP() * (*m);
                          return m;
                      }),
                      allow_raw_pointers())
            .function("vectorD", optional_override([](const Eigen::LDLT<MatrixXf>& self) {
                          return new MatrixXf(self.vectorD());
                      }),
                      allow_raw_pointers())
            .function("solve", optional_override([](const Eigen::LDLT<MatrixXf>& self, MatrixXf& b) {
                          return new MatrixXf(self.solve(b));
                      }),
                      allow_raw_pointers());

    class_<Eigen::LLT<MatrixXf>>("LLTf")
            .function("matrixL", optional_override([](const Eigen::LLT<MatrixXf>& self) {
                          return new MatrixXf(self.matrixL());
                      }),
                      allow_raw_pointers())
            .function("solve", optional_override([](const Eigen::LLT<MatrixXf>& self, MatrixXf& b) {
                          return new MatrixXf(self.solve(b));
                      }),
                      allow_raw_pointers());

    class_<FullPivLUXf>("FullPivLUf")
            .function("solve", optional_override([](FullPivLUXf& self, MatrixXf& b) {
                          MatrixXf* m = new MatrixXf(self.solve(b));
                          return m;
                      }),
                      allow_raw_pointers())
            .function("matrixLU", optional_override([](FullPivLUXf& self) {
                          return new MatrixXf(self.matrixLU());
                      }),
                      allow_raw_pointers())
            .function("rank", &FullPivLUXf::rank)
            .function("isInvertible", &FullPivLUXf::isInvertible)
            .function("isSurjective", &FullPivLUXf::isSurjective)
            .function("isInjective", &FullPivLUXf::isInjective)
            .function("determinant", &FullPivLUXf::determinant)
            .function("rcond", &FullPivLUXf::rcond)
            .function("dimensionOfKernel", &FullPivLUXf::dimensionOfKernel)
            .function("inverse", optional_override([](const FullPivLUXf& self) {
                          return new MatrixXf(self.inverse());
                      }),
                      allow_raw_pointers())
            .function("matrixL", optional_override([](const FullPivLUXf& self) {
                          int rows = self.rows();
                          int cols = self.cols();
                          MatrixXf* L = new MatrixXf(rows, rows);
                          L->setIdentity();
                          L->block(0, 0, rows, cols).triangularView<Eigen::StrictlyLower>() = self.matrixLU();
                          return L;
                      }),
                      allow_raw_pointers())
            .function("matrixU", optional_override([](const FullPivLUXf& self) {
                          MatrixXf* U = new MatrixXf(self.matrixLU().triangularView<Eigen::Upper>());
                          return U;
                      }),
                      allow_raw_pointers())
            .function("permutationQ", optional_override([](const FullPivLUXf& self) {
                          return new MatrixXf(self.permutationQ());
                      }),
                      allow_raw_pointers())
            .function("permutationP", optional_override([](const FullPivLUXf& self) {
                          return new MatrixXf(self.permutationP());
                      }),
                      allow_raw_pointers());

    class_<PartialPivLUXf>("PartialPivLUf")
            .function("solve", optional_override([](PartialPivLUXf& self, MatrixXf& b) {
                          return new MatrixXf(self.solve(b));
                      }),
                      allow_raw_pointers())
            .function("matrixLU", optional_override([](PartialPivLUXf& self) {
                          return new MatrixXf(self.matrixLU());
                      }),
                      allow_raw_pointers())
            .function("determinant", &PartialPivLUXf::determinant)
            .function("rcond", &PartialPivLUXf::rcond)
            .function("inverse", optional_override([](const PartialPivLUXf& self) {
                          return new MatrixXf(self.inverse());
                      }),
                      allow_raw_pointers())
            .function("matrixL", optional_override([](const PartialPivLUXf& self) {
                          int rows = self.rows();
                          int cols = self.cols();
                          MatrixXf* L = new MatrixXf(rows, rows);
                          L->setIdentity();
                          L->block(0, 0, rows, cols).triangularView<Eigen::StrictlyLower>() = self.matrixLU();
                          return L;
                      }),
                      allow_raw_pointers())
            .function("matrixU", optional_override([](const PartialPivLUXf& self) {
                          MatrixXf* U = new MatrixXf(self.matrixLU().triangularView<Eigen::Upper>());
                          return U;
                      }),
                      allow_raw_pointers())
            .function("permutationP", optional_override([](const PartialPivLUXf& self) {
                          return new MatrixXf(self.permutationP());
                      }),
                      allow_raw_pointers());

    class_<HouseholderQRXf>("HouseholderQRf")
            .function("solve", optional_override([](HouseholderQRXf& self, MatrixXf& b) {
                          return new MatrixXf(self.solve(b));
                      }),
                      allow_raw_pointers())
            .function("householderQ", optional_override([](HouseholderQRXf& self) {
                          return new MatrixXf(self.householderQ());
                      }),
                      allow_raw_pointers());

    class_<ColPivHouseholderQRXf>("ColPivHouseholderQRf")
            .function("solve", optional_override([](ColPivHouseholderQRXf& self, MatrixXf& b) {
                          return new MatrixXf(self.solve(b));
                      }),
                      allow_raw_pointers())
            .function("colsPermutation", optional_override([](const ColPivHouseholderQRXf& self) {
                          return new MatrixXf(self.colsPermutation());
                      }),
                      allow_raw_pointers())
            .function("householderQ", optional_override([](ColPivHouseholderQRXf& self) {
                          return new MatrixXf(self.householderQ());
                      }),
                      allow_raw_pointers())
            .function("matrixR", optional_override([](ColPivHouseholderQRXf& self) {
                          return new MatrixXf(self.matrixR());
                      }),
                      allow_raw_pointers())
            .function("matrixQR", optional_override([](ColPivHouseholderQRXf& self) {
                          return new MatrixXf(self.matrixQR());
                      }),
                      allow_raw_pointers())
            .function("rank", &ColPivHouseholderQRXf::rank)
            .function("isInvertible", &ColPivHouseholderQRXf::isInvertible)
            .function("isSurjective", &ColPivHouseholderQRXf::isSurjective)
            .function("isInjective", &ColPivHouseholderQRXf::isInjective)
            .function("logAbsDeterminant", &ColPivHouseholderQRXf::logAbsDeterminant)
            .function("dimensionOfKernel", &ColPivHouseholderQRXf::dimensionOfKernel)
            .function("inverse", optional_override([](const ColPivHouseholderQRXf& self) {
                          return new MatrixXf(self.inverse());
                      }),
                      allow_raw_pointers());

    class_<FullPivHouseholderQRXf>("FullPivHouseholderQRf")
            .function("solve", optional_override([](FullPivHouseholderQRXf& self, MatrixXf& b) {
                          return new MatrixXf(self.solve(b));
                      }),
                      allow_raw_pointers())
            .function("colsPermutation", optional_override([](const FullPivHouseholderQRXf& self) {
                          return new MatrixXf(self.colsPermutation());
                      }),
                      allow_raw_pointers())
            .function("matrixQ", optional_override([](FullPivHouseholderQRXf& self) {
                          return new MatrixXf(self.matrixQ());
                      }),
                      allow_raw_pointers())
            .function("matrixQR", optional_override([](FullPivHouseholderQRXf& self) {
                          return new MatrixXf(self.matrixQR());
                      }),
                      allow_raw_pointers())
            .function("rank", &FullPivHouseholderQRXf::rank)
            .function("isInvertible", &FullPivHouseholderQRXf::isInvertible)
            .function("isSurjective", &FullPivHouseholderQRXf::isSurjective)
            .function("isInjective", &FullPivHouseholderQRXf::isInjective)
            .function("logAbsDeterminant", &FullPivHouseholderQRXf::logAbsDeterminant)
            .function("dimensionOfKernel", &FullPivHouseholderQRXf::dimensionOfKernel)
            .function("inverse", optional_override([](const FullPivHouseholderQRXf& self) {
                          return new MatrixXf(self.inverse());
                      }),
                      allow_raw_pointers());

    class_<BDCSVDXf>("BDCSVDf")
            .function("solve", optional_override([](BDCSVDXf& self, MatrixXf& b) {
                          return new MatrixXf(self.solve(b));
                      }),
                      allow_raw_pointers())
            .function("matrixU", optional_override([](BDCSVDXf& self) {
                          return new MatrixXf(self.matrixU());
                      }),
                      allow_raw_pointers())
            .function("matrixV", optional_override([](BDCSVDXf& self) {
                          return new MatrixXf(self.matrixV());
                      }),
                      allow_raw_pointers())
            .function("nonzeroSingularValues", &BDCSVDXf::nonzeroSingularValues)
            .function("rank", &BDCSVDXf::rank)
            .function("singularValues", optional_override([](const BDCSVDXf& self) {
                          return new Eigen::MatrixXf(self.singularValues());
                      }),
                      allow_raw_pointers());

    class_<SelfAdjointEigenSolverXf>("SelfAdjointEigenSolverf")
            .function("eigenvalues", optional_override([](SelfAdjointEigenSolverXf& self) {
                          return new MatrixXf(self.eigenvalues());
                      }),
                      allow_raw_pointers())
            .function("eigenvectors", optional_override([](SelfAdjointEigenSolverXf& self) {
                          return new MatrixXf(self.eigenvectors());
                      }),
                      allow_raw_pointers());
#endif// DECOMPOSITIONS
}