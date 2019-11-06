import lldb


# Synthetic view for static_vector objects
class static_vector_SynthChildProvider:
  def __init__( self, valobj, dict ):
    self.valobj = valobj

  def num_children( self ):
    size = self.valobj.GetChildMemberWithName( '_size' ).GetValueAsUnsigned()
    return size + 2

  def has_children( self ):
    return self.num_children() != 0

  def get_child_at_index( self, index ):
    if index < 0:
      return None
    if index >= self.num_children():
      return None
    if index == self.num_children() - 2:
      return self.valobj.GetChildMemberWithName( '_arr' )
    if index == self.num_children() - 1:
      return self.valobj.GetChildMemberWithName( '_size' )
    arr = self.valobj.GetChildMemberWithName( '_arr' )
    return arr.GetChildAtIndex(index)

  def update( self ):
    pass

