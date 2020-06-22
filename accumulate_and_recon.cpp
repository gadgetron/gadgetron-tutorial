#include <gadgetron/Node.h>
#include <gadgetron/mri_core_data.h>
#include <gadgetron/hoNDFFT.h>
#include <gadgetron/mri_core_utility.h>
#include <gadgetron/mri_core_coil_map_estimation.h>

using namespace Gadgetron;

class AccumulateAndReconGadget : public Core::ChannelGadget<Core::Acquisition> {

    public: 
        AccumulateAndReconGadget(const Core::Context& context, const Core::GadgetProperties& props) : Core::ChannelGadget<Core::Acquisition>(context, props){

        }


    virtual void process(Core::InputChannel<Core::Acquisition>& in, Core::OutputChannel& out){

        auto matrix_size = this->header.encoding.front().encodedSpace.matrixSize;
        hoNDArray<std::complex<float>> data(matrix_size.x,matrix_size.y,matrix_size.z, this->header.acquisitionSystemInformation->receiverChannels);


        for (auto [acq_header, data, trajectory] : in) {
            using namespace Gadgetron::Indexing;

            //data(:,acq_header.idx.kspace_encode_step_1,acq_header.idx.kspace_encode_step_2,:) = data;
            data(slice,acq_header.idx.kspace_encode_step_1,acq_header.idx.kspace_encode_step_2,slice) = data;

            if (acq_header.isFlagSet(ISMRMRD::ISMRMRD_ACQ_LAST_IN_ENCODE_STEP2)){

                auto image_data = FFT::ifft3c(data);
                auto coil_map = coil_map_Inati(data);
                image_data = coil_combine(data,coil_map,4);

                auto image_header = image_header_from_acquisition(acq_header,this->header,image_data);
                out.push(image_header,image_data);

                data.fill(0.0f);

            }

        }
    }

    
};

GADGETRON_GADGET_EXPORT(AccumulateAndReconGadget);