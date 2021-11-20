// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once
#include <stdio.h>

class mcParticle;

class mcPhaseSpaceIO
{
public:
	// Types
	enum open_mode_e { NOT_OPEN, READING, WRITING };

	// Ctor
	mcPhaseSpaceIO();
	// Dtor
	virtual ~mcPhaseSpaceIO();

	// Accessors
	unsigned int numParticles()const { return preamble_.n_particles; }
	unsigned int latch()const { return latch_; }
	bool bitReg(unsigned int)const;
	unsigned int birthReg()const;

	// Modifiers
	void open(const char *filename, open_mode_e open_mode);
	void close();
	mcParticle read();
	void write(const mcParticle& particle);

private:
	// Types
	struct record_short_t {
		unsigned int latch_;
		float energy_;
		float x_, y_;
		float u_, v_;
		float weight_times_sign_w_; //negative sign of w is possible
	};
	struct record_long_t {
		unsigned int latch_;
		float energy_;
		float x_, y_;
		float u_, v_;
		float weight_times_sign_w_; //negative sign of w is possible
		float zlast_;
	};
	struct preamble_t {
		//char mode[5];
		unsigned int n_particles;
		unsigned int n_photons;
		float ekin_max;
		float e_min_electrons;
		float n_incident_particles;
	};
	enum mode_e { UNKNOWN, MODE0, MODE2 };

	void stripZLast(); //shorten the record size from record_long_t to record_short_t
	void reset();
	void flush();  // write next particle
	void fill();   // read next particle
	unsigned startData();

	// Data
	// Constants
	static const unsigned int PHOTON_LATCH;
	static const unsigned int NEGATRON_LATCH;
	static const unsigned int POSITRON_LATCH;
	static const unsigned int PROTON_LATCH;
	static const unsigned int NEUTRON_LATCH;
	static const unsigned int BUFFER_SIZE;

	// R/W Buffers
	record_short_t *short_buffer;
	record_long_t *long_buffer;

	// Data used for all operations
	open_mode_e open_mode_;
	mode_e mode_;

	unsigned int buf_pos;  // Number of records in buffer
	unsigned int read_pos; // Position to read next in the buffer
	FILE *file;
	preamble_t preamble_;

	// Keep latch of the last read particle
	unsigned int latch_;
};
