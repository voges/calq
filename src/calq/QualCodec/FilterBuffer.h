class FilterBuffer {
private:
    std::vector<double> kernel;
    std::vector<double> buffer;
    size_t bufferPos; //Index pointer to oldest value

public:
    //New activity score in pipeline
    void push (double activityScore);

    //Calculate filter score at offset position
    double filter const ();

    //Initialize buffer and 
    FilterBuffer(const std::function<double, size_t, size_t>& kernelBuilder, size_t kernelSize);

    //Create dummy buffer
    FilterBuffer();

    //Buffer size
    size_t getSize() const;

    //DIstance between buffer center and borders
    size_t getOffset() const;
}
