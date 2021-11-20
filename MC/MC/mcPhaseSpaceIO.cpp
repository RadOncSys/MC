#include "mcPhaseSpaceIO.h"
#include "mcParticle.h"
#include "mcDefs.h"
#include <float.h>

const unsigned int mcPhaseSpaceIO::PHOTON_LATCH = 0x00000000;
const unsigned int mcPhaseSpaceIO::NEGATRON_LATCH = 0x40000000;
const unsigned int mcPhaseSpaceIO::POSITRON_LATCH = 0x20000000;
const unsigned int mcPhaseSpaceIO::PROTON_LATCH = 0x80000000;
const unsigned int mcPhaseSpaceIO::NEUTRON_LATCH = 0x40000000;
const unsigned int mcPhaseSpaceIO::BUFFER_SIZE = 1024;

mcPhaseSpaceIO::mcPhaseSpaceIO()
	:open_mode_(NOT_OPEN)
	, mode_(UNKNOWN)
	, buf_pos(0)
	, read_pos(0)
	, file(nullptr)
{
	short_buffer = new record_short_t[BUFFER_SIZE];
	long_buffer = new record_long_t[BUFFER_SIZE];
	reset();
}

mcPhaseSpaceIO::~mcPhaseSpaceIO()
{
	close();
	delete[] short_buffer;
	delete[] long_buffer;
}

void mcPhaseSpaceIO::reset()
{
	mode_ = UNKNOWN;
	if (file) {
		fclose(file);
		file = nullptr;
	}
	open_mode_ = NOT_OPEN;
	buf_pos = 0;
	read_pos = 0;
	preamble_.n_particles = 0;
	preamble_.n_photons = 0;
	preamble_.ekin_max = (float)0.;
	preamble_.e_min_electrons = FLT_MAX;
	preamble_.n_incident_particles = float(preamble_.n_particles);
}

void mcPhaseSpaceIO::close()
{
	// If we were writing, flush the buffers and update the preamble
	if (open_mode_ == WRITING) {
		flush();
		fseek(file, 5, SEEK_SET);
		fwrite(&preamble_, sizeof(preamble_t), 1, file);
	}
	reset();
}

void mcPhaseSpaceIO::stripZLast()
{
	unsigned int i;
	for (i = 0; i < buf_pos; i++) {
		short_buffer[i] = *reinterpret_cast<record_short_t *>(long_buffer + i);
	}
}

bool mcPhaseSpaceIO::bitReg(unsigned int i)const
{
	unsigned int mask = 0x00000001;
	mask <<= i;
	return (latch_ & mask) > 0;
}

unsigned int mcPhaseSpaceIO::birthReg()const
{
	unsigned int reg = latch_ & 0x1F000000;
	reg >>= 24;
	return reg;
}

void mcPhaseSpaceIO::flush()
{
	if (open_mode_ == WRITING) {
		fwrite(short_buffer, sizeof(record_short_t), buf_pos, file); // ONLY MODE0
		buf_pos = 0;
	}
}

void mcPhaseSpaceIO::fill()
{
	if (open_mode_ == READING) {

		// Try to read the BUFFER_SIZE record, and set buf_size to the actual number
		// of records read
		if (mode_ == MODE2) {
			buf_pos = (unsigned)fread(long_buffer, sizeof(record_long_t), BUFFER_SIZE, file);
			stripZLast();
		}
		else {
			buf_pos = (unsigned)fread(short_buffer, sizeof(record_short_t), BUFFER_SIZE, file);
		}
		read_pos = 0;
		// If we've reached the EOF, rewind the stream
		if (feof(file)) {
			fseek(file, startData(), SEEK_SET);
		}
	}
}

unsigned mcPhaseSpaceIO::startData()
{
	unsigned offset = (mode_ == MODE0) ? 3 : 7;
	return offset + 5 + sizeof(preamble_t);
}

void mcPhaseSpaceIO::open(const char *filename, open_mode_e open_mode)
{
	const char *om;
	switch (open_mode) {
	case READING:
		om = "rb";
		break;
	case WRITING:
		om = "wb";
		break;
	default:
		throw std::exception("Unknown opening mode");
	}
	reset();

	if (fopen_s(&file, filename, om) != 0)
		throw std::exception("Cannot open file");
	open_mode_ = open_mode;
	if (open_mode_ == WRITING) mode_ = MODE0; // Can write only MODE0 files
	else { // Read preamble
		char modekey[5];
		fread(modekey, sizeof(char), 5, file);
		fread(&preamble_, sizeof(preamble_t), 1, file);

		if (!strncmp(modekey, "MODE0", 5)) mode_ = MODE0;
		else if (!strncmp(modekey, "MODE2", 5)) mode_ = MODE2;
		char dummy[7];
		switch (mode_) {
		case MODE0:
			fread(dummy, 3, 1, file);
			break;
		case MODE2:
			fread(dummy, 7, 1, file);
			break;
		default:
			throw std::exception("Format of phase space file is different from MODE0 or MODE2");
		}
	}
}

mcParticle mcPhaseSpaceIO::read()
{
	if (open_mode_ != READING) throw std::exception("File is not open for reading");
	if (read_pos >= buf_pos) fill();
	if (!buf_pos) fill();
	if (!buf_pos) throw std::exception("No records in file");

	const record_short_t& rec = short_buffer[read_pos];
	read_pos++;

	latch_ = rec.latch_;
	geomVector3D r(rec.x_, rec.y_, 0.);
	double uz2 = 1. - rec.u_*rec.u_ - rec.v_*rec.v_;
	if (uz2 < 0) uz2 = 0;
	if (uz2 < 0.) throw std::exception("Bad direction");
	double uz = sqrt(uz2);
	if (rec.weight_times_sign_w_ < 0.) uz = -uz;
	geomVector3D u(rec.u_, rec.v_, uz);

	mc_particle_t type = MCP_PHOTON;
	double ke = rec.energy_;
	int pq = 0;
	double weight = fabs(rec.weight_times_sign_w_);

	if (rec.latch_ & NEGATRON_LATCH)
	{
		type = MCP_NEGATRON;
		ke -= EMASS;
		pq = -1;
	}
	else if (rec.latch_ & POSITRON_LATCH)
	{
		type = MCP_POSITRON;
		ke -= EMASS;
		pq = 1;
	}
	else if (rec.latch_ & PROTON_LATCH)
	{
		type = MCP_PROTON;
		ke -= PMASS;
		pq = 1;
	}
	else if (rec.latch_ & NEUTRON_LATCH)
	{
		type = MCP_NEUTRON;
	}
	mcParticle particle(type, pq, ke, r, u);
	particle.weight = weight;
	return particle;
}

void mcPhaseSpaceIO::write(const mcParticle& p)
{
	if (open_mode_ != WRITING) throw std::exception("File is not open for writing");
	if (!preamble_.n_particles) {
		fwrite("MODE0", sizeof(char), 5, file);
		fwrite(&preamble_, sizeof(preamble_t), 1, file);
		char dummy[7];
		fwrite(dummy, 3, 1, file);
	}
	if (buf_pos >= BUFFER_SIZE) {
		flush();
	}
	preamble_.n_particles++;
	preamble_.n_incident_particles = float(preamble_.n_particles);
	if (preamble_.ekin_max < p.ke) preamble_.ekin_max = (float)p.ke;
	record_short_t rec;
	rec.x_ = (float)p.p.x();
	rec.y_ = (float)p.p.y();
	rec.u_ = (float)p.u.x();
	rec.v_ = (float)p.u.y();
	rec.weight_times_sign_w_ = (float)p.weight;
	if (p.u.z() < 0.) rec.weight_times_sign_w_ = -rec.weight_times_sign_w_;

	switch (p.t)
	{
	case MCP_PHOTON:
		rec.latch_ = PHOTON_LATCH;
		rec.energy_ = (float)p.ke;
		preamble_.n_photons += 1;
		break;
	case MCP_POSITRON:
		rec.latch_ = POSITRON_LATCH;
		rec.energy_ = float(p.ke + EMASS);
		if (preamble_.e_min_electrons > float(p.ke + TWICE_EMASS)) preamble_.e_min_electrons = float(p.ke + TWICE_EMASS);
		break;
	case MCP_NEGATRON:
		rec.latch_ = NEGATRON_LATCH;
		rec.energy_ = float(p.ke + EMASS);
		if (preamble_.e_min_electrons > float(p.ke)) preamble_.e_min_electrons = float(p.ke);
		break;
	case MCP_PROTON:
		rec.latch_ = PROTON_LATCH;
		rec.energy_ = float(p.ke + PMASS);
		//if(preamble_.e_min_electrons > p.ke) preamble_.e_min_electrons = p.ke;
		break;
	case MCP_NEUTRON:
		rec.latch_ = NEUTRON_LATCH;
		rec.energy_ = float(p.ke);
		//if(preamble_.e_min_electrons > p.ke) preamble_.e_min_electrons = p.ke;
		break;
	default:
		throw std::exception("Unknown particle type");
	}
	short_buffer[buf_pos] = rec;
	buf_pos++;
}
