namespace safetics
{
    class Motion
    {
    public:
        Motion();
        ~Motion();
        int MouseControl();
        void LoadDHParameters();
        void setTCPMatrix();
        void ForwardKinematics();
        void InverseKinematics();
    };

    void create2DArray(int** arr, int row, int col);
    void TMatrix(double** T, double a, double alpha, double d, double theta);
}