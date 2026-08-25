function mesh = make_ert_mesh_fvm_rect(varargin)
%MAKE_ERT_MESH_FVM_RECT Convenience wrapper for refined rectilinear FVM meshes.

    mesh = make_ert_mesh_hex8(varargin{:});
    mesh.meshStyle = "rect_refined";
end
