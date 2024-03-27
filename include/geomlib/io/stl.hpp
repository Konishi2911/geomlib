#pragma once
#include <stdexcept>
#ifndef GEOMLIB_IO_STL_HPP
#define GEOMLIB_IO_STL_HPP

#include "geomlib/surface/polygon_mesh.hpp"
#include "geomlib/io/binary_utility.hpp"

namespace geomlib {

struct StlError {
    struct IlegalPolygon: std::runtime_error {
        IlegalPolygon(): std::runtime_error("STL cannot be constructed from the mesh containing non-triangular faces.") {}
    };
};

/// @brief	A class representing a STL (Standard Triangle Language) model.
/// @remark This class provides only export as binary file. 
///			A feature to load a .stl file and export as ascii file will be added in the future.
struct Stl {
public:
	Stl(PolygonMesh&& faces, std::string&& header = "");

	template<typename Os> 
    void write_binary(Os& os);

private:
	const PolygonMesh _mesh;
	const std::string _header;

    static auto __resize_header(std::string&& str) noexcept -> std::string;
};


inline Stl::Stl(PolygonMesh&& mesh, std::string&& header):
    _mesh(std::move(mesh)), 
    _header(__resize_header(std::move(header))) 
{
    auto&& faces = this->_mesh.subfaces();
    bool is_triangle = true;
    for (auto i = 0u; i < this->_mesh.n_faces(); ++i) {
        is_triangle &= faces[0].verts().size() == 3;
    }

    if (!is_triangle) {
        throw StlError::IlegalPolygon();
    } 
}


template<typename Os> 
void Stl::write_binary(Os& os) {
	if (os) {
		auto buf = std::vector<uint8_t>();
		buf.reserve(84 + this->_mesh.n_faces() * (3 * 4 * 8));

		// Header section
		{
			auto header_bin = io::as_binary(this->_header);
			buf.insert(buf.end(), header_bin.begin(), header_bin.end());
		}

		// Face count section
		{ 
			auto n_faces = static_cast<uint32_t>(this->_mesh.n_faces());
			auto n_faces_bin = io::as_binary(n_faces);
			buf.insert(buf.end(), n_faces_bin.begin(), n_faces_bin.end());
		}

		// Face section
		for (const auto& face: this->_mesh.subfaces()) {
			// normal vector
			{
				auto n = face.normal();
				auto n_float = std::array<float, 3>();
				std::transform(
					n.cbegin(), n.cend(), n_float.begin(), 
					[](const double& d){ return static_cast<float>(d); }
				);

				auto n_bin = io::into_binary(n_float.cbegin(), n_float.cend());
				buf.insert(buf.end(), n_bin.begin(), n_bin.end());
			}

			// vertices
			for (auto&& vert: face.verts()) {
				auto vert_float = std::array<float, 3>();
				std::transform(
					vert.get().cbegin(), vert.get().cend(), vert_float.begin(), 
					[](const double& d){ return static_cast<float>(d); }
				);

				auto vert_bin = io::into_binary(vert_float.cbegin(), vert_float.cend());
				buf.insert(buf.end(), vert_bin.begin(), vert_bin.end());
			}

			// attribute byte count
			buf.insert(buf.end(), {0, 0});
		}
		os.write(reinterpret_cast<char*>(buf.data()), buf.size());
	}
}

inline auto Stl::__resize_header(std::string&& str) noexcept -> std::string { 
	str.resize(80);
	return std::move(str);
}

}

#endif